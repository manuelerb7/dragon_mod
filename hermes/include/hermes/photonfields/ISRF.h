#ifndef HERMES_ISRF_H
#define HERMES_ISRF_H

#include <array>

#include "hermes/photonfields/PhotonField.h"
#include "hermes/FITSWrapper.h"

namespace hermes { namespace photonfields {
/**
 * \addtogroup PhotonFields
 * @{
 */

class ISRF : public PhotonField {
  private:
	const static int freqR1 = 200;
	const static int freqR2 = 680;
	const static int freqR3 = 331;
	
	std::unique_ptr<FITSFile> ffile;
	int n_r, n_z, n_wvl;
	double min_r, min_z, min_wvl;
	double delta_r, delta_z, delta_wvl;
	//std::array<double, freqR1 + freqR2 + freqR3> logwavelenghts;
	//std::array<double, 30> r_id = {
	//    0.0,  0.2,  0.5,  1.0,  1.5,  2.0,  2.5,  3.0, 3.5, 4.0, 4.5,
	//    5.0,  5.5,  6.0,  6.5,  7.0,  7.5,  8.0,  8.5, 9.0, 9.5, 10.0,
	//    11.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0};  // in kpc (30)
	//std::array<double, 24> z_id = {0.0,  0.1,  0.2, 0.3, 0.4,  0.5,  0.6,
	//                               0.8,  1.0,  1.2, 1.5, 2.0,  2.5,  3.0,
	//                               4.0,  5.0,  6.0, 8.0, 10.0, 12.0, 15.0,
	//                               20.0, 25.0, 30.0};  // in kpc (24)
	std::vector<double> logwavelenghts;
	std::vector<double> r_id;  // in kpc (30)
	std::vector<double> z_id;  // in kpc (24)
	std::vector<float> isrf;

	void buildEnergyRange();

	double getISRF(std::size_t ir, std::size_t iz, std::size_t ifreq) const;
	double getFitsISRF(std::size_t ir, std::size_t iz, std::size_t ifreq) const;

	void loadFrequencyAxis();
	void loadISRF();
	void loadFitsISRF(const std::string &isrf_filename);

  public:
	ISRF(const int isrf_mode);
	
	int read_mode;
	
	std::size_t getSize() const;
	QEnergyDensity getEnergyDensity(const QLength &r, const QLength &z,
	                                const QEnergy &E_photon) const;
	QEnergyDensity getEnergyDensity(const Vector3QLength &pos,
	                                const QEnergy &E_photon) const override;
	QEnergyDensity getEnergyDensity(const Vector3QLength &pos_,
	                                std::size_t iE_) const override;
};

/** @}*/
}}  // namespace hermes::photonfields

#endif  // HERMES_ISRF_H
