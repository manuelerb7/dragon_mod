#ifdef HERMES_HAVE_CFITSIO

#include "hermes/photonfields/ISRF.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

#include "hermes/Common.h"

namespace hermes { namespace photonfields {

std::string str(const int &n) {
	std::stringstream ss;
	ss << std::setfill('0') << std::setw(3) << n;
	return ss.str();
}

ISRF::ISRF(const int isrf_mode) {
	
	read_mode = isrf_mode;
	
	auto logWavelenghtToFrequency = [](double lambda) {
		return c_light / (std::pow(10, lambda) * micrometre);
	};
	

	// Spares steps
	// setEnergyScaleFactor(1.1); // 145 steps
	

	// Alternative (slow), all available energy steps
	/*
	double scaling = std::pow(static_cast<double>(
	            getEndEnergy()/getStartEnergy()),
	        1.0/logwavelenghts.size()); // ~1.01
	setEnergyScaleFactor(scaling); // 1211 steps
	*/

	
	
	if (read_mode == 0) {
		r_id = {
		    0.0,  0.2,  0.5,  1.0,  1.5,  2.0,  2.5,  3.0, 3.5, 4.0, 4.5,
		    5.0,  5.5,  6.0,  6.5,  7.0,  7.5,  8.0,  8.5, 9.0, 9.5, 10.0,
		    11.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0};  // in kpc (30)
		
		z_id = {0.0,  0.1,  0.2, 0.3, 0.4,  0.5,  0.6,
		    0.8,  1.0,  1.2, 1.5, 2.0,  2.5,  3.0,
		    4.0,  5.0,  6.0, 8.0, 10.0, 12.0, 15.0,
		    20.0, 25.0, 30.0};
		
		loadFrequencyAxis();
		
		loadISRF();
	}
	else if (read_mode == 1) {
	
		std::string isrf_filename = "MilkyWay_DR0.5_DZ0.1_DPHI10_RMAX20_ZMAX5_galprop_format.fits.gz";
		
		loadFitsISRF(isrf_filename);
	}
	
	else {
		std::stringstream ss;
		ss << "hermes: error: only ISRF mode 0 and 1 available";
		throw std::runtime_error(ss.str());
	}
	
	setStartEnergy(logWavelenghtToFrequency(logwavelenghts.back()) * h_planck);
	setEndEnergy(logWavelenghtToFrequency(logwavelenghts.front()) * h_planck);
	setEnergyScaleFactor(1.05);
			
	buildEnergyRange();
	
}

void ISRF::buildEnergyRange() {
	const double scaling = getEnergyScaleFactor();
	const QEnergy E_start = getStartEnergy();
	const QEnergy E_end = getEndEnergy();

	for (QEnergy E = E_start; E < E_end; E = E * scaling)
		energyRange.push_back(E);
}

void ISRF::loadFrequencyAxis() {
	double logwl = log10(0.01);  // micron
	for (size_t i = 0; i < freqR1; ++i) {
		logwavelenghts.push_back(logwl);
		logwl += 0.01;
	}
	for (size_t i = 0; i < freqR2; ++i) {
		logwavelenghts.push_back(logwl);
		logwl += 0.0025;
	}
	for (size_t i = 0; i < freqR3; ++i) {
		logwavelenghts.push_back(logwl);
		logwl += 0.01;
	}
}

std::size_t ISRF::getSize() const { return isrf.size(); }

void ISRF::loadISRF() {
	const int max_num_of_char_in_a_line = 512;
	const int num_of_header_lines = 1;

	int n = 0;
	for (auto i : r_id) {
		for (auto j : z_id) {
			std::ostringstream name;
			name << "RadiationField/Vernetto16/spectrum_r"
			     << str(static_cast<int>(i * 10)) << "_z"
			     << str(static_cast<int>(j * 10)) << ".dat";
			std::string filename = getDataPath(name.str());

			std::ifstream fin(filename.c_str());
			if (!fin) {
				std::stringstream ss;
				ss << "hermes: error: File " << filename << " not found";
				throw std::runtime_error(ss.str());
			}
			for (std::size_t k = 0; k < num_of_header_lines; ++k) {
				fin.ignore(max_num_of_char_in_a_line, '\n');
			}
			while (!fin.eof()) {
				float f_, e_;
				fin >> f_ >> e_;
				if (!fin.eof()) isrf.push_back(e_);
			}
			n++;
		}
	}
	assert(isrf.size() == r_id.size() * z_id.size() * logwavelenghts.size());
}

void ISRF::loadFitsISRF(const std::string &isrf_filename) {
	
	ffile = std::make_unique<FITSFile>(FITSFile(getDataPath("RadiationField/" + isrf_filename)));
	ffile->openFile(FITS::READ);

	n_r = ffile->readKeyValueAsInt("NAXIS1");
	n_z = ffile->readKeyValueAsInt("NAXIS2");
	n_wvl = ffile->readKeyValueAsInt("NAXIS3");

	min_r = ffile->readKeyValueAsDouble("CRVAL1");
	delta_r = ffile->readKeyValueAsDouble("CDELT1");
	min_z = ffile->readKeyValueAsDouble("CRVAL2");
	delta_z = ffile->readKeyValueAsDouble("CDELT2");
	min_wvl = ffile->readKeyValueAsDouble("CRVAL3");
	delta_wvl = ffile->readKeyValueAsDouble("CDELT3");
	
	for (size_t i = 0; i < n_r; ++i) {
		r_id.push_back(min_r + i*delta_r);
	}
	
	for (size_t i = 0; i < n_z; ++i) {
		z_id.push_back(min_z + i*delta_z);
	}
	
	for (size_t i = 0; i < n_wvl; ++i) {
		logwavelenghts.push_back(min_wvl + i*delta_wvl);
	}
	
	//setStartEnergy(logWavelenghtToFrequency(logwavelenghts.back()) * h_planck);
	//setEndEnergy(logWavelenghtToFrequency(logwavelenghts.front()) * h_planck);
	//setEnergyScaleFactor(1.05);
		
	//buildEnergyRange();

	int firstElement = 1;
	int nElements = n_r * n_z * n_wvl;
	isrf = ffile->readImageAsFloat(firstElement, nElements);
	
	assert(isrf.size() == r_id.size() * z_id.size() * logwavelenghts.size());
}


double ISRF::getISRF(std::size_t ir, std::size_t iz, std::size_t imu) const {
	std::size_t i = imu + iz * logwavelenghts.size() +
	                ir * (logwavelenghts.size() * z_id.size());
	return isrf[i];
}

double ISRF::getFitsISRF(std::size_t ir, std::size_t iz, std::size_t imu) const {
	std::size_t i = ir + iz * r_id.size() + imu * r_id.size() * z_id.size();
	return isrf[i];
}

QEnergyDensity ISRF::getEnergyDensity(const Vector3QLength &pos,
                                      std::size_t iE) const {
	QLength r = sqrt(pos.x * pos.x + pos.y * pos.y);
	QLength z = pos.z;
	QEnergy E = energyRange[iE];
	// TODO(adundovi): not implemented
	return getEnergyDensity(r, z, E);
}

QEnergyDensity ISRF::getEnergyDensity(const Vector3QLength &pos,
                                      const QEnergy &E_photon) const {
	QLength r = sqrt(pos.x * pos.x + pos.y * pos.y);
	QLength z = pos.z;
	return getEnergyDensity(r, z, E_photon);
}

QEnergyDensity ISRF::getEnergyDensity(const QLength &r, const QLength &z,
                                      const QEnergy &E_photon) const {
	double r_ = static_cast<double>(r / 1_kpc);
	double z_ = static_cast<double>(fabs(z) / 1_kpc);
	double f_mu =
	    static_cast<double>(h_planck * c_light / E_photon / (micrometre));
	double logf_ = std::log10(f_mu);

	if (r_ < r_id.front() || r_ > r_id.back()) return 0;
	if (z_ < z_id.front() || z_ > z_id.back()) return 0;
	if (logf_ < logwavelenghts.front() || logf_ > logwavelenghts.back())
		return 0;

	std::size_t ir =
	    std::lower_bound(r_id.begin(), r_id.end(), r_) - r_id.begin();
	std::size_t iz =
	    std::lower_bound(z_id.begin(), z_id.end(), z_) - z_id.begin();
	std::size_t ifreq =
	    std::lower_bound(logwavelenghts.begin(), logwavelenghts.end(), logf_) -
	    logwavelenghts.begin() - 1;

	if (ir == r_id.size()) return 0;
	if (iz == z_id.size()) return 0;
	if (ifreq == logwavelenghts.size()) return 0;

	double r_d = (r_ - r_id[ir]) / (r_id[ir + 1] - r_id[ir]);
	double z_d = (z_ - z_id[iz]) / (z_id[iz + 1] - z_id[iz]);
	double f_d = (logf_ - logwavelenghts[ifreq]) /
	             (logwavelenghts[ifreq + 1] - logwavelenghts[ifreq]);

	/*
	if (!(r_d >= 0 && r_d <= 1))
	    return 0;
	if (!(z_d >= 0 && z_d <= 1))
	    return 0;
	if (!(f_d >= 0 && f_d <= 1))
	    return 0;
	*/
	
	double c_00, c_01, c_10, c_11;
	
	if (read_mode == 0) {
		c_00 = getISRF(ir, iz, ifreq) * (1. - r_d) + getISRF(ir + 1, iz, ifreq) * r_d;
		c_01 = getISRF(ir, iz, ifreq + 1) * (1. - r_d) + getISRF(ir + 1, iz, ifreq + 1) * r_d;
		c_10 = getISRF(ir, iz + 1, ifreq) * (1. - r_d) + getISRF(ir + 1, iz + 1, ifreq) * r_d;
		c_11 = getISRF(ir, iz + 1, ifreq + 1) * (1. - r_d) + getISRF(ir + 1, iz + 1, ifreq + 1) * r_d;
	} else if (read_mode == 1) {
		c_00 = getFitsISRF(ir, iz, ifreq) * (1. - r_d) + getFitsISRF(ir + 1, iz, ifreq) * r_d;
		c_01 = getFitsISRF(ir, iz, ifreq + 1) * (1. - r_d) + getFitsISRF(ir + 1, iz, ifreq + 1) * r_d;
		c_10 = getFitsISRF(ir, iz + 1, ifreq) * (1. - r_d) + getFitsISRF(ir + 1, iz + 1, ifreq) * r_d;
		c_11 = getFitsISRF(ir, iz + 1, ifreq + 1) * (1. - r_d) + getFitsISRF(ir + 1, iz + 1, ifreq + 1) * r_d;
	}
	double c_0 = c_00 * (1. - z_d) + c_10 * z_d;
	double c_1 = c_01 * (1. - z_d) + c_11 * z_d;

	double c = c_0 * (1. - f_d) + c_1 * f_d;

	return c * 1_eV / 1_cm3;
}

}}  // namespace hermes::photonfields

#endif  // HERMES_HAVE_CFITSIO
