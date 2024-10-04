// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::EclHysteresisTwoPhaseLawParams
 */
#ifndef OPM_ECL_HYSTERESIS_TWO_PHASE_LAW_PARAMS_HPP
#define OPM_ECL_HYSTERESIS_TWO_PHASE_LAW_PARAMS_HPP

#include <opm/input/eclipse/EclipseState/WagHysteresisConfig.hpp>

#include <opm/material/common/EnsureFinalized.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsConfig.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>
#include <opm/material/fluidmatrixinteractions/EclHysteresisConfig.hpp>

#include <cassert>
#include <cmath>
#include <memory>
namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief A default implementation of the parameters for the material law which
 *        implements the ECL relative permeability and capillary pressure hysteresis
 */
template <class EffLawT>
class EclHysteresisTwoPhaseLawParams : public EnsureFinalized
{
    using EffLawParams = typename EffLawT::Params;
    using Scalar = typename EffLawParams::Traits::Scalar;

public:
    using Traits = typename EffLawParams::Traits;

    static EclHysteresisTwoPhaseLawParams serializationTestObject()
    {
        EclHysteresisTwoPhaseLawParams<EffLawT> result;
        result.dynamic_.deltaSwImbKrn_ = 1.0;
        //result.deltaSwImbKrw_ = 1.0;
        result.dynamic_.Sncrt_ = 2.0;
        result.dynamic_.Swcrt_ = 2.5;
        result.dynamic_.initialImb_ = true;
        result.dynamic_.pcSwMic_ = 3.0;
        result.dynamic_.krnSwMdc_ = 4.0;
        result.dynamic_.krwSwMdc_ = 4.5;
        result.dynamic_.KrndHy_ = 5.0;
        result.dynamic_.KrwdHy_ = 6.0;

        return result;
    }

    /*!
     * \brief Calculate all dependent quantities once the independent
     *        quantities of the parameter object have been set.
     */
    void finalize()
    {
        if (config().enableHysteresis()) {
            if (config().krHysteresisModel() == 2 || config().krHysteresisModel() == 3 || config().krHysteresisModel() == 4 || config().pcHysteresisModel() == 0) {
                C_ = 1.0/(Sncri_ - Sncrd_ + 1.0e-12) - 1.0/(Snmaxd_ - Sncrd_);
                curvatureCapPrs_ =  config().curvatureCapPrs();
            }
            if (config().krHysteresisModel() == 4) {
                Cw_ = 1.0/(Swcri_ - Swcrd_ + 1.0e-12) - 1.0/(Swmaxd_ - Swcrd_);
                Krwd_sncri_ = EffLawT::twoPhaseSatKrw(drainageParams(), 1 - Sncri());
            }
            updateDynamicParams_();
        }
        EnsureFinalized :: finalize();
    }

    /*!
     * \brief Set the hysteresis configuration object.
     */
    void setConfig(std::shared_ptr<EclHysteresisConfig> value)
    { config_ = *value; }

    /*!
     * \brief Returns the hysteresis configuration object.
     */
    const EclHysteresisConfig& config() const
    { return config_; }

    /*!
     * \brief Set the WAG-hysteresis configuration object.
     */
    void setWagConfig(std::shared_ptr<WagHysteresisConfig::WagHysteresisConfigRecord> value)
    {
        wagConfig_ = value;
        dynamic_.cTransf_ = wagConfig().wagLandsParam();
    }

    /*!
     * \brief Returns the WAG-hysteresis configuration object.
     */
    const WagHysteresisConfig::WagHysteresisConfigRecord& wagConfig() const
    { return *wagConfig_; }

    /*!
     * \brief Sets the parameters used for the drainage curve
     */
    void setDrainageParams(const EffLawParams& value,
                           const EclEpsScalingPointsInfo<Scalar>& info,
                           EclTwoPhaseSystemType twoPhaseSystem)
    {
        drainageParams_ = value;

        oilWaterSystem_ = (twoPhaseSystem == EclTwoPhaseSystemType::OilWater);
        gasOilSystem_ = (twoPhaseSystem == EclTwoPhaseSystemType::GasOil);

        if (!config().enableHysteresis())
            return;
        if (config().krHysteresisModel() == 2 || config().krHysteresisModel() == 3 || config().krHysteresisModel() == 4 || config().pcHysteresisModel() == 0 || gasOilHysteresisWAG()) {
            Swco_ = info.Swl;
            if (twoPhaseSystem == EclTwoPhaseSystemType::GasOil) {
                Sncrd_ = info.Sgcr + info.Swl;
                Swcrd_ = info.Sogcr;
                Snmaxd_ = info.Sgu + info.Swl;
                KrndMax_ = EffLawT::twoPhaseSatKrn(drainageParams(), 1.0-Snmaxd_);
            }
            else if (twoPhaseSystem == EclTwoPhaseSystemType::GasWater) {
                Sncrd_ = info.Sgcr;
                Swcrd_ = info.Swcr;
                Snmaxd_ = info.Sgu;
                KrndMax_ = EffLawT::twoPhaseSatKrn(drainageParams(), 1.0-Snmaxd_);
            }
            else {
                assert(twoPhaseSystem == EclTwoPhaseSystemType::OilWater);
                Sncrd_ = info.Sowcr;
                Swcrd_ = info.Swcr;
                Snmaxd_ = 1.0 - info.Swl - info.Sgl;
                KrndMax_ = EffLawT::twoPhaseSatKrn(drainageParams(), 1.0-Snmaxd_);
            }
        }

        if (config().krHysteresisModel() == 4) {
            //Swco_ = info.Swl;
            if (twoPhaseSystem == EclTwoPhaseSystemType::GasOil) {
                Swmaxd_ = 1.0 - info.Sgl - info.Swl; 
                KrwdMax_ = EffLawT::twoPhaseSatKrw(drainageParams(), Swmaxd_);
            }
            else if (twoPhaseSystem == EclTwoPhaseSystemType::GasWater) {
                Swmaxd_ = info.Swu;
                KrwdMax_ = EffLawT::twoPhaseSatKrw(drainageParams(), Swmaxd_);
            }
            else {
                assert(twoPhaseSystem == EclTwoPhaseSystemType::OilWater);
                Swmaxd_ = info.Swu;
                KrwdMax_ = EffLawT::twoPhaseSatKrw(drainageParams(), Swmaxd_);
            }
        }

        // Additional Killough hysteresis model for pc
        if (config().pcHysteresisModel() == 0) {
            if (twoPhaseSystem == EclTwoPhaseSystemType::GasOil) {
                pcmaxd_ = info.maxPcgo;
            } else if (twoPhaseSystem == EclTwoPhaseSystemType::GasWater) {
                pcmaxd_ = info.maxPcgo + info.maxPcow;
            }
            else {
                assert(twoPhaseSystem == EclTwoPhaseSystemType::OilWater);
                pcmaxd_ = -17.0; // At this point 'info.maxPcow' holds pre-swatinit value ...;
            }
        }

        // For WAG hysteresis, assume initial state along primary drainage curve.
        if (gasOilHysteresisWAG()) {
            dynamic_.swatImbStart_ = Swco_;
            dynamic_.swatImbStartNxt_ = -1.0; // Trigger check for saturation gt Swco at first update ...
            dynamic_.cTransf_ = wagConfig().wagLandsParam();
            dynamic_.krnSwDrainStart_ = Sncrd_;
            dynamic_.krnSwDrainStartNxt_ = Sncrd_;
            dynamic_.krnImbStart_ = 0.0;
            dynamic_.krnImbStartNxt_ = 0.0;
            dynamic_.krnDrainStart_ = 0.0;
            dynamic_.krnDrainStartNxt_ = 0.0;
            dynamic_.isDrain_ = true;
            dynamic_.wasDrain_ = true;
            dynamic_.krnSwImbStart_ = Sncrd_;
            dynamic_.SncrtWAG_ = Sncrd_;
            dynamic_.nState_ = 1;
        }
    }

    /*!
     * \brief Returns the parameters used for the drainage curve
     */
    const EffLawParams& drainageParams() const
    { return drainageParams_; }

    EffLawParams& drainageParams()
    { return drainageParams_; }

    /*!
     * \brief Sets the parameters used for the imbibition curve
     */
    void setImbibitionParams(const EffLawParams& value,
                             const EclEpsScalingPointsInfo<Scalar>& info,
                             EclTwoPhaseSystemType twoPhaseSystem)
    {
        imbibitionParams_ = value;

        if (!config().enableHysteresis())
            return;

        // Store critical nonwetting saturation
        if (twoPhaseSystem == EclTwoPhaseSystemType::GasOil) {
            Sncri_ = info.Sgcr + info.Swl;
            Swcri_ = info.Sogcr;
        }
        else if (twoPhaseSystem == EclTwoPhaseSystemType::GasWater) {
            Sncri_ = info.Sgcr;
            Swcri_ = info.Swcr;
        }
        else {
            assert(twoPhaseSystem == EclTwoPhaseSystemType::OilWater);
            Sncri_ = info.Sowcr;
            Swcri_ = info.Swcr;
        }

        // Killough hysteresis model for pc
        if (config().pcHysteresisModel() == 0) {
            if (twoPhaseSystem == EclTwoPhaseSystemType::GasOil) {
                Swmaxi_ = 1.0 - info.Sgl - info.Swl;
                pcmaxi_ = info.maxPcgo;
            } else if (twoPhaseSystem == EclTwoPhaseSystemType::GasWater) {
                Swmaxi_ = 1.0 - info.Sgl;
                pcmaxi_ = info.maxPcgo + info.maxPcow;
            }
            else {
                assert(twoPhaseSystem == EclTwoPhaseSystemType::OilWater);
                Swmaxi_ = info.Swu;
                pcmaxi_ = info.maxPcow;
            }
        }

        if (config().krHysteresisModel() == 4) {
            Krwi_snmax_ = EffLawT::twoPhaseSatKrw(imbibitionParams(), 1 - Snmaxd());
            Krwi_snrmax_ = EffLawT::twoPhaseSatKrw(imbibitionParams(), 1 - Sncri());
        }
    }

    /*!
     * \brief Returns the parameters used for the imbibition curve
     */
    const EffLawParams& imbibitionParams() const
    { return imbibitionParams_; }

    EffLawParams& imbibitionParams()
    { return imbibitionParams_; }

    /*!
     * \brief Get the saturation of the wetting phase where the last switch from the main
     *        drainage curve to imbibition happend on the capillary pressure curve.
     */
    Scalar pcSwMdc() const
    { return dynamic_.pcSwMdc_; }

    Scalar pcSwMic() const
    { return dynamic_.pcSwMic_; }

    /*!
     * \brief Status of initial process.
     */
    bool initialImb() const
    { return dynamic_.initialImb_; }

    /*!
     * \brief Set the saturation of the wetting phase where the last switch from the main
     *        drainage curve (MDC) to imbibition happend on the relperm curve for the
     *        wetting phase.
     */
    void setKrwSwMdc(Scalar value)
    { dynamic_.krwSwMdc_ = value; };

    /*!
     * \brief Get the saturation of the wetting phase where the last switch from the main
     *        drainage curve to imbibition happend on the relperm curve for the
     *        wetting phase.
     */
    Scalar krwSwMdc() const
    { return dynamic_.krwSwMdc_; };

    /*!
     * \brief Set the saturation of the wetting phase where the last switch from the main
     *        drainage curve (MDC) to imbibition happend on the relperm curve for the
     *        non-wetting phase.
     */
    void setKrnSwMdc(Scalar value)
    { dynamic_.krnSwMdc_ = value; }

    /*!
     * \brief Get the saturation of the wetting phase where the last switch from the main
     *        drainage curve to imbibition happend on the relperm curve for the
     *        non-wetting phase.
     */
    Scalar krnSwMdc() const
    { return dynamic_.krnSwMdc_; }

    /*!
     * \brief Sets the saturation value which must be added if krw is calculated using
     *        the imbibition curve.
     *
     * This means that krw(Sw) = krw_drainage(Sw) if Sw < SwMdc and
     * krw(Sw) = krw_imbibition(Sw + Sw_shift,krw) else
     */
    //void setDeltaSwImbKrw(Scalar value)
    //{ deltaSwImbKrw_ = value; }

    /*!
     * \brief Returns the saturation value which must be added if krw is calculated using
     *        the imbibition curve.
     *
     * This means that krw(Sw) = krw_drainage(Sw) if Sw < SwMdc and
     * krw(Sw) = krw_imbibition(Sw + Sw_shift,krw) else
     */
    //Scalar deltaSwImbKrw() const
    //{ return deltaSwImbKrw_; }

    /*!
     * \brief Sets the saturation value which must be added if krn is calculated using
     *        the imbibition curve.
     *
     * This means that krn(Sw) = krn_drainage(Sw) if Sw < SwMdc and
     * krn(Sw) = krn_imbibition(Sw + Sw_shift,krn) else
     */
    void setDeltaSwImbKrn(Scalar value)
    { dynamic_.deltaSwImbKrn_ = value; }

    /*!
     * \brief Returns the saturation value which must be added if krn is calculated using
     *        the imbibition curve.
     *
     * This means that krn(Sw) = krn_drainage(Sw) if Sw < SwMdc and
     * krn(Sw) = krn_imbibition(Sw + Sw_shift,krn) else
     */
    Scalar deltaSwImbKrn() const
    { return dynamic_.deltaSwImbKrn_; }


    Scalar Swcri() const
    { return Swcri_; }

    Scalar Swcrd() const
    { return Swcrd_; }

    Scalar Swmaxi() const
    { return Swmaxi_; }

    Scalar Sncri() const
    { return Sncri_; }

    Scalar Sncrd() const
    { return Sncrd_; }

    Scalar Sncrt() const
    { return dynamic_.Sncrt_; }

    Scalar Swcrt() const
    { return dynamic_.Swcrt_; }

    Scalar SnTrapped(bool maximumTrapping) const
    {
        if(!maximumTrapping && dynamic_.isDrain_)
            return 0.0;

        // For Killough the trapped saturation is already computed
        if( config().krHysteresisModel() > 1 )
            return dynamic_.Sncrt_;
        else // For Carlson we use the shift to compute it from the critial saturation
            return Sncri_ + dynamic_.deltaSwImbKrn_;
    }

    Scalar SnStranded(Scalar sg, Scalar krg) const {
        const Scalar sn = EffLawT::twoPhaseSatKrnInv(drainageParams_, krg);
        return sg - (1.0 - sn) + Sncrd_;
    }

    Scalar SwTrapped() const
    {
        if( config().krHysteresisModel() == 0 || config().krHysteresisModel() == 2)
            return Swcrd_;
        
        if( config().krHysteresisModel() == 1 || config().krHysteresisModel() == 3)
            return Swcri_;
        
        // For Killough the trapped saturation is already computed
        if( config().krHysteresisModel() == 4 )
            return dynamic_.Swcrt_;
        
        return 0.0;
        //else // For Carlson we use the shift to compute it from the critial saturation
        //    return Swcri_ + deltaSwImbKrw_;
    }

    Scalar SncrtWAG() const
    { return dynamic_.SncrtWAG_; }

    Scalar Snmaxd() const
    { return Snmaxd_; }

    Scalar Swmaxd() const
    { return Swmaxd_; }

    Scalar Snhy() const
    { return 1.0 - dynamic_.krnSwMdc_; }

    Scalar Swhy() const
    { return dynamic_.krwSwMdc_; }

    Scalar Swco() const
    { return Swco_; }

    Scalar krnWght() const
    { return dynamic_.KrndHy_/KrndMax_; }

    Scalar krwWght() const
    {
        // a = 1 (deltaKrw)^a Formulation according to KILLOUGH 1976
        Scalar deltaKrw = Krwi_snrmax() - Krwd_sncri();
        Scalar Krwi_snr = Krwd_sncrt() + deltaKrw * (Sncrt() / max(1e-12, Sncri()));
        return (Krwi_snr - KrwdHy()) / ( Krwi_snrmax() - Krwi_snmax());
    }

    Scalar krwdMax() const
    { 
        return KrwdMax_; }

    Scalar KrwdHy() const
    {
        return dynamic_.KrwdHy_;
    }


    Scalar Krwd_sncri() const
    {
        return Krwd_sncri_;
    }

    Scalar Krwi_snmax() const
    {
        return Krwi_snmax_;
    }

    Scalar Krwi_snrmax() const
    {
        return Krwi_snrmax_;
    }

    Scalar Krwd_sncrt() const
    {
        return dynamic_.Krwd_sncrt_;
    }

    Scalar pcWght() const // Aligning pci and pcd at Swir
    {
        if (pcmaxd_ < 0.0)
            return EffLawT::twoPhaseSatPcnw(drainageParams(), 0.0)/(pcmaxi_+1e-6);
        else
            return pcmaxd_/(pcmaxi_+1e-6);
    }

    Scalar curvatureCapPrs() const
    { return curvatureCapPrs_;}

    bool gasOilHysteresisWAG() const
    { return (config().enableWagHysteresis() && gasOilSystem_ && wagConfig().wagGasFlag()) ; }

    Scalar reductionDrain() const
    { return std::pow(Swco_/(dynamic_.swatImbStart_+tolWAG_*wagConfig().wagWaterThresholdSaturation()), wagConfig().wagSecondaryDrainageReduction());}

    Scalar reductionDrainNxt() const
    { return std::pow(Swco_/(dynamic_.swatImbStartNxt_+tolWAG_*wagConfig().wagWaterThresholdSaturation()), wagConfig().wagSecondaryDrainageReduction());}

    bool threePhaseState() const
    { return (dynamic_.swatImbStart_ > (Swco_ + wagConfig().wagWaterThresholdSaturation()) ); }

    Scalar nState() const
    { return dynamic_.nState_;}

    Scalar krnSwDrainRevert() const
    { return dynamic_.krnSwDrainRevert_;}

    Scalar krnDrainStart() const
    { return dynamic_.krnDrainStart_;}

    Scalar krnDrainStartNxt() const
    { return dynamic_.krnDrainStartNxt_;}

    Scalar krnImbStart() const
    { return dynamic_.krnImbStart_;}

    Scalar krnImbStartNxt() const
    { return dynamic_.krnImbStartNxt_;}

    Scalar krnSwWAG() const
    { return dynamic_.krnSwWAG_;}

    Scalar krnSwDrainStart() const
    { return dynamic_.krnSwDrainStart_;}

    Scalar krnSwDrainStartNxt() const
    { return dynamic_.krnSwDrainStartNxt_;}

    Scalar krnSwImbStart() const
    { return dynamic_.krnSwImbStart_;}

    Scalar tolWAG() const
    { return tolWAG_;}

    template <class Evaluation>
    Evaluation computeSwf(const Evaluation& Sw)  const
    {
        Evaluation SgT = 1.0 - Sw - SncrtWAG(); // Sg-Sg_crit_trapped
        Scalar SgCut = wagConfig().wagImbCurveLinearFraction()*(Snhy()- SncrtWAG());
        Evaluation Swf = 1.0;
        //Scalar C = wagConfig().wagLandsParam();
        Scalar C = dynamic_.cTransf_;

        if (SgT > SgCut) {
            Swf -= (Sncrd() + 0.5*( SgT + Opm::sqrt( SgT*SgT + 4.0/C*SgT))); // 1-Sgf
        }
        else {
            SgCut = std::max(Scalar(0.000001), SgCut);
            Scalar SgCutValue = 0.5*( SgCut + Opm::sqrt( SgCut*SgCut + 4.0/C*SgCut));
            Scalar SgCutSlope = SgCutValue/SgCut;
            SgT *= SgCutSlope;
            Swf -= (Sncrd() + SgT);
        }

        return Swf;
    }

    template <class Evaluation>
    Evaluation computeKrImbWAG(const Evaluation& Sw)  const
    {
        Evaluation Swf = Sw;
        if (dynamic_.nState_ <= 2)  // Skipping for "higher order" curves seems consistent with benchmark, further investigations needed ...
            Swf = computeSwf(Sw);
        if (Swf <= dynamic_.krnSwDrainStart_) { // Use secondary drainage curve
            Evaluation Krg = EffLawT::twoPhaseSatKrn(drainageParams_, Swf);
            Evaluation KrgImb2 = (Krg-dynamic_.krnDrainStart_)*reductionDrain() + dynamic_.krnImbStart_;
            return KrgImb2;
        }
        else { // Fallback to primary drainage curve
            Evaluation Sn = Sncrd_;
            if (Swf < 1.0-dynamic_.SncrtWAG_) {
                // Notation: Sn.. = Sg.. + Swco
                Evaluation dd = (1.0-dynamic_.krnSwImbStart_ - Sncrd_) / (1.0-dynamic_.krnSwDrainStart_ - dynamic_.SncrtWAG_);
                Sn += (1.0-Swf-dynamic_.SncrtWAG_)*dd;
            }
            Evaluation KrgDrn1 = EffLawT::twoPhaseSatKrn(drainageParams_, 1.0 - Sn);
            return KrgDrn1;
        }
    }

    /*!
     * \brief Notify the hysteresis law that a given wetting-phase saturation has been seen
     *
     * This updates the scanning curves and the imbibition<->drainage reversal points as
     * appropriate.
     */
    bool update(Scalar pcSw, Scalar krwSw, Scalar krnSw)
    {
        bool updateParams = false;

        if (config().pcHysteresisModel() == 0 && pcSw < dynamic_.pcSwMdc_) {
            if (dynamic_.pcSwMdc_ == 2.0 && pcSw+1.0e-6 < Swcrd_ && oilWaterSystem_) {
               dynamic_.initialImb_ = true;
            }
            dynamic_.pcSwMdc_ = pcSw;
            updateParams = true;
        }

        if (dynamic_.initialImb_ && pcSw > dynamic_.pcSwMic_) {
            dynamic_.pcSwMic_ = pcSw;
            updateParams = true;
        }

        if (krnSw < dynamic_.krnSwMdc_) {
            dynamic_.krnSwMdc_ = krnSw;
            dynamic_.KrndHy_ = EffLawT::twoPhaseSatKrn(drainageParams(), dynamic_.krnSwMdc_);
            if (config().krHysteresisModel() == 4) {
                dynamic_.KrwdHy_ = EffLawT::twoPhaseSatKrw(drainageParams(), dynamic_.krnSwMdc_);
            }
            updateParams = true;
        }
        if (krwSw > dynamic_.krwSwMdc_) {
            dynamic_.krwSwMdc_ = krwSw; // Only used for output at the moment
        }

        // for non WAG hysteresis we still keep track of the process
        // for output purpose.
        if (!gasOilHysteresisWAG()) {
            this->dynamic_.isDrain_ = (krnSw <= this->dynamic_.krnSwMdc_);
        } else {
            dynamic_.wasDrain_ = dynamic_.isDrain_;

            if (dynamic_.swatImbStartNxt_ < 0.0) { // Initial check ...
                dynamic_.swatImbStartNxt_ = std::max(Swco_, Swco_ + krnSw - krwSw);
                // check if we are in threephase state sw > swco + tolWag and so > tolWag
                // (sw = swco + krnSw - krwSw and so = krwSw for oil/gas params)
                if ( (dynamic_.swatImbStartNxt_ > Swco_ + tolWAG_) && krwSw > tolWAG_) {
                    dynamic_.swatImbStart_ = dynamic_.swatImbStartNxt_;
                    dynamic_.krnSwWAG_ = krnSw;
                    dynamic_.krnSwDrainStartNxt_ = dynamic_.krnSwWAG_;
                    dynamic_.krnSwDrainStart_ = dynamic_.krnSwDrainStartNxt_;
                    dynamic_.wasDrain_ = false; // Signal start from threephase state ...
                }
            }

            if (dynamic_.isDrain_) {
                if (krnSw <= dynamic_.krnSwWAG_+tolWAG_) { // continue along drainage curve
                    dynamic_.krnSwWAG_ = std::min(krnSw, dynamic_.krnSwWAG_);
                    dynamic_.krnSwDrainRevert_ = dynamic_.krnSwWAG_;
                    updateParams = true;
                }
                else { // start new imbibition curve
                    dynamic_.isDrain_ = false;
                    dynamic_.krnSwWAG_ = krnSw;
                    updateParams = true;
                }
            }
            else {
                if (krnSw >= dynamic_.krnSwWAG_-tolWAG_) { // continue along imbibition curve
                    dynamic_.krnSwWAG_ = std::max(krnSw, dynamic_.krnSwWAG_);
                    dynamic_.krnSwDrainStartNxt_ = dynamic_.krnSwWAG_;
                    dynamic_.swatImbStartNxt_ = std::max(dynamic_.swatImbStartNxt_, Swco_ + krnSw - krwSw);
                    updateParams = true;
                }
                else { // start new drainage curve
                    dynamic_.isDrain_ = true;
                    dynamic_.krnSwDrainStart_ = dynamic_.krnSwDrainStartNxt_;
                    dynamic_.swatImbStart_ = dynamic_.swatImbStartNxt_;
                    dynamic_.krnSwWAG_ = krnSw;
                    updateParams = true;
                }
            }

        }

        if (updateParams)
            updateDynamicParams_();

        return updateParams;
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        // only serializes dynamic state - see update() and updateDynamic_()
        serializer(dynamic_.deltaSwImbKrn_);
        //serializer(deltaSwImbKrw_);
        serializer(dynamic_.Sncrt_);
        serializer(dynamic_.Swcrt_);
        serializer(dynamic_.initialImb_);
        serializer(dynamic_.pcSwMic_);
        serializer(dynamic_.krnSwMdc_);
        serializer(dynamic_.krwSwMdc_);
        serializer(dynamic_.KrndHy_);
        serializer(dynamic_.KrwdHy_);
    }

    bool operator==(const EclHysteresisTwoPhaseLawParams& rhs) const
    {
        return this->dynamic_.deltaSwImbKrn_ == rhs.dynamic_.deltaSwImbKrn_ &&
               //this->deltaSwImbKrw_ == rhs.deltaSwImbKrw_ &&
               this->dynamic_.Sncrt_ == rhs.dynamic_.Sncrt_ &&
               this->dynamic_.Swcrt_ == rhs.dynamic_.Swcrt_ &&
               this->dynamic_.initialImb_ == rhs.dynamic_.initialImb_ &&
               this->dynamic_.pcSwMic_ == rhs.dynamic_.pcSwMic_ &&
               this->dynamic_.krnSwMdc_ == rhs.dynamic_.krnSwMdc_ &&
               this->dynamic_.krwSwMdc_ == rhs.dynamic_.krwSwMdc_ &&
               this->dynamic_.KrndHy_ == rhs.dynamic_.KrndHy_ &&
               this->dynamic_.KrwdHy_ == rhs.dynamic_.KrwdHy_;
    }

private:
    void updateDynamicParams_() const
    {
        // calculate the saturation deltas for the relative permeabilities
        //if (false) { // we dont support Carlson for wetting phase hysteresis
            //Scalar krwMdcDrainage = EffLawT::twoPhaseSatKrw(drainageParams(), krwSwMdc_);
            //Scalar SwKrwMdcImbibition = EffLawT::twoPhaseSatKrwInv(imbibitionParams(), krwMdcDrainage);
            //deltaSwImbKrw_ = SwKrwMdcImbibition - krwSwMdc_;
        //}   

        if (config().krHysteresisModel() == 0 || config().krHysteresisModel() == 1) {
            Scalar krnMdcDrainage = EffLawT::twoPhaseSatKrn(drainageParams(), dynamic_.krnSwMdc_);
            Scalar SwKrnMdcImbibition = EffLawT::twoPhaseSatKrnInv(imbibitionParams(), krnMdcDrainage);
            dynamic_.deltaSwImbKrn_ = SwKrnMdcImbibition - dynamic_.krnSwMdc_;
        }

        // Scalar pcMdcDrainage = EffLawT::twoPhaseSatPcnw(drainageParams(), pcSwMdc_);
        // Scalar SwPcMdcImbibition = EffLawT::twoPhaseSatPcnwInv(imbibitionParams(), pcMdcDrainage);
        // deltaSwImbPc_ = SwPcMdcImbibition - pcSwMdc_;

        if (config().krHysteresisModel() == 2 || config().krHysteresisModel() == 3 || config().krHysteresisModel() == 4 || config().pcHysteresisModel() == 0) {
            Scalar Snhy = 1.0 - dynamic_.krnSwMdc_;
            if (Snhy > Sncrd_) {
                dynamic_.Sncrt_ = Sncrd_ + (Snhy - Sncrd_)/((1.0+config().modParamTrapped()*(Snmaxd_-Snhy)) + C_*(Snhy - Sncrd_));
            }
            else
            {
                dynamic_.Sncrt_ = Sncrd_;
            }
        }

        if (config().krHysteresisModel() == 4) {
            Scalar Swhy = dynamic_.krnSwMdc_;
            if (Swhy >= Swcrd_) {
                dynamic_.Swcrt_ = Swcrd_ + (Swhy - Swcrd_)/((1.0+config().modParamTrapped()*(Swmaxd_-Swhy)) + Cw_*(Swhy - Swcrd_));
            } else
            {
                dynamic_.Swcrt_ = Swcrd_;
            }
            dynamic_.Krwd_sncrt_ = EffLawT::twoPhaseSatKrw(drainageParams(), 1 - Sncrt());
        }


        if (gasOilHysteresisWAG()) {
            if (dynamic_.isDrain_ && dynamic_.krnSwMdc_ == dynamic_.krnSwWAG_) {
                Scalar Snhy = 1.0 - dynamic_.krnSwMdc_;
                dynamic_.SncrtWAG_ = Sncrd_;
                if (Snhy > Sncrd_) {
                    dynamic_.SncrtWAG_ += (Snhy - Sncrd_)/(1.0+config().modParamTrapped()*(Snmaxd_-Snhy) + wagConfig().wagLandsParam()*(Snhy - Sncrd_));
                }
            }

            if (dynamic_.isDrain_ && (1.0-dynamic_.krnSwDrainRevert_) > dynamic_.SncrtWAG_) { //Reversal from drain to imb
                dynamic_.cTransf_ = 1.0/(dynamic_.SncrtWAG_-Sncrd_ + 1.0e-12) - 1.0/(1.0-dynamic_.krnSwDrainRevert_-Sncrd_);
            }

            if (!dynamic_.wasDrain_ && dynamic_.isDrain_) { // Start of new drainage cycle
                if (threePhaseState() || dynamic_.nState_>1) { // Never return to primary (two-phase) state after leaving
                    dynamic_.nState_ += 1;
                    dynamic_.krnDrainStart_ = EffLawT::twoPhaseSatKrn(drainageParams(), dynamic_.krnSwDrainStart_);
                    dynamic_.krnImbStart_ = dynamic_.krnImbStartNxt_;
                    // Scanning shift for primary drainage
                    dynamic_.krnSwImbStart_ = EffLawT::twoPhaseSatKrnInv(drainageParams(), dynamic_.krnImbStart_);
                }
            }

            if (!dynamic_.wasDrain_ && !dynamic_.isDrain_) { //Moving along current imb curve
                dynamic_.krnDrainStartNxt_ = EffLawT::twoPhaseSatKrn(drainageParams(), dynamic_.krnSwWAG_);
                if (threePhaseState()) {
                    dynamic_.krnImbStartNxt_ = computeKrImbWAG(dynamic_.krnSwWAG_);
                }
                else {
                    Scalar swf = computeSwf(dynamic_.krnSwWAG_);
                    dynamic_.krnImbStartNxt_ = EffLawT::twoPhaseSatKrn(drainageParams(), swf);
                }
            }

        }

    }

    struct DynamicHysteresisState
    {
        // offsets added to wetting phase saturation uf using the imbibition curves need to
        // be used to calculate the wetting phase relperm, the non-wetting phase relperm and
        // the capillary pressure
        //Scalar deltaSwImbKrw_{};
        Scalar deltaSwImbKrn_{};
        //Scalar deltaSwImbPc_;

        Scalar Swcrt_{}; // trapped wetting phase saturation
        Scalar Krwd_sncrt_{};
        Scalar Sncrt_{}; // trapped non-wetting phase saturation
        Scalar SncrtWAG_{};
        Scalar cTransf_{}; // Modified Lands constant used for free gas calculations to obtain consistent scanning curve
                           //  when reversion to imb occurs above historical maximum gas saturation (i.e. Sw > krwSwMdc_).
        int nState_{};                 // Number of cycles. Primary cycle is nState_=1.
        Scalar krnDrainStart_{};       // Primary (input) relperm evaluated at start of current drainage curve.
        Scalar krnDrainStartNxt_{};    // Primary (input) relperm evaluated at start of next drainage curve.
        Scalar krnImbStart_{};         // Relperm at start of current drainage curve (end of previous imb curve).
        Scalar krnImbStartNxt_{};      // Relperm at start of next drainage curve (end of current imb curve).
        Scalar krnSwImbStart_{};       // Saturation value where primary drainage relperm equals krnImbStart_

        // largest wettinging phase saturation which is on the main-drainage curve. These are
        // three different values because the sourounding code can choose to use different
        // definitions for the saturations for different quantities
        Scalar krwSwMdc_{-2.0};
        Scalar krnSwMdc_{2.0};
        Scalar pcSwMdc_{2.0};

        // largest wettinging phase saturation along main imbibition curve
        Scalar pcSwMic_{1.0};
        // Initial process is imbibition (for initial saturations at or below critical drainage saturation)
        bool initialImb_{false};

        Scalar KrndHy_{};  // Krn_drain(1-krnSwMdc_)
        Scalar KrwdHy_{};

        bool isDrain_{true};           // Status is either drainage or imbibition
        bool wasDrain_{};              // Previous status.

        Scalar swatImbStart_{};        // Water saturation at start of current drainage curve (end of previous imb curve).
        Scalar swatImbStartNxt_{};     // Water saturation at start of next drainage curve (end of current imb curve).

        Scalar krnSwWAG_{2.0};         // Saturation value after latest completed timestep.
        Scalar krnSwDrainStart_{-2.0}; // Saturation value at start of current drainage curve (end of previous imb curve).
        Scalar krnSwDrainStartNxt_{};  // Saturation value at start of current drainage curve (end of previous imb curve).
        Scalar krnSwDrainRevert_{2.0}; // Saturation value at end of current drainage curve.

    };

    EclHysteresisConfig config_{};
    std::shared_ptr<WagHysteresisConfig::WagHysteresisConfigRecord> wagConfig_{};
    EffLawParams imbibitionParams_{};
    EffLawParams drainageParams_{};



    bool oilWaterSystem_{false};
    bool gasOilSystem_{false};

    // the following uses the conventions of the Eclipse technical description:
    //
    // Sncrd_: critical non-wetting phase saturation for the drainage curve
    // Sncri_: critical non-wetting phase saturation for the imbibition curve
    // Swcri_: critical wetting phase saturation for the imbibition curve
    // Swcrd_: critical wetting phase saturation for the drainage curve
    // Swmaxi_; maximum wetting phase saturation for the imbibition curve
    // Snmaxd_: non-wetting phase saturation where the non-wetting relperm reaches its
    //          maximum on the drainage curve
    // C_: factor required to calculate the trapped non-wetting phase saturation using
    //     the Killough approach
    // Cw_: factor required to calculate the trapped wetting phase saturation using
    //     the Killough approach
    // Swcod_: connate water saturation value used for wag hysteresis (2. drainage)
    Scalar Sncrd_{};
    Scalar Sncri_{};
    Scalar Swcri_{};
    Scalar Swcrd_{};
    Scalar Swmaxi_{};
    Scalar Snmaxd_{};
    Scalar Swmaxd_{};
    Scalar C_{};

    Scalar KrndMax_{}; // Krn_drain(Snmaxd_)
    Scalar KrwdMax_{}; // Krw_drain(Swmaxd_)

    // For wetting hysterese Killough
    Scalar Cw_{};
    Scalar Krwd_sncri_{};
    Scalar Krwi_snmax_{};
    Scalar Krwi_snrmax_{};

    Scalar pcmaxd_{};  // max pc for drain
    Scalar pcmaxi_{};  // max pc for imb

    Scalar curvatureCapPrs_{}; // curvature parameter used for capillary pressure hysteresis

    // Used for WAG hysteresis
    Scalar Swco_{};                // Connate water.

    Scalar tolWAG_{0.001};

    mutable DynamicHysteresisState dynamic_;
};

} // namespace Opm

#endif
