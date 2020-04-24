#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

namespace trajectoryStateTransform {

  using namespace SurfaceSideDefinition;

  PTrajectoryStateOnDet persistentState(const TrajectoryStateOnSurface& ts, unsigned int detid) {
    int surfaceSide = static_cast<int>(ts.surfaceSide());
    auto pt = ts.globalMomentum().perp();

    if (ts.hasError()) {
      AlgebraicSymMatrix55 const& m = ts.localError().matrix();

      int dim = 5;  /// should check if corresponds to m
      float localErrors[15];

      int k = 0;
      for (int i = 0; i < dim; i++) {
        for (int j = 0; j <= i; j++) {
          localErrors[k++] = m(i, j);
        }
      }
      return PTrajectoryStateOnDet(ts.localParameters(), pt, localErrors, detid, surfaceSide);
    }
    return PTrajectoryStateOnDet(ts.localParameters(), pt, detid, surfaceSide);
  }

  TrajectoryStateOnSurface transientState(const PTrajectoryStateOnDet& ts,
                                          const Surface* surface,
                                          const MagneticField* field) {
    AlgebraicSymMatrix55 m;
    bool errInv = true;
    if (ts.hasError()) {
      errInv = false;
      int dim = 5;
      int k = 0;
      for (int i = 0; i < dim; i++) {
        for (int j = 0; j <= i; j++) {
          m(i, j) = ts.error(k++);  // NOTE: here we do a cast float => double.
        }
      }
    }

    return TrajectoryStateOnSurface(ts.parameters(),
                                    errInv ? LocalTrajectoryError(InvalidError()) : LocalTrajectoryError(m),
                                    *surface,
                                    field,
                                    static_cast<SurfaceSide>(ts.surfaceSide()));
  }

  FreeTrajectoryState initialFreeState(const reco::Track& tk, const MagneticField* field, bool withErr) {
    Basic3DVector<float> pos(tk.vertex());
    GlobalPoint gpos(pos);
    Basic3DVector<float> mom(tk.momentum());
    GlobalVector gmom(mom);
    GlobalTrajectoryParameters par(gpos, gmom, tk.charge(), field);
    if (!withErr)
      return FreeTrajectoryState(par);
    CurvilinearTrajectoryError err(tk.covariance());
    return FreeTrajectoryState(par, err);
  }

  FreeTrajectoryState initialFreeStateTTrack(const TTTrack< Ref_Phase2TrackerDigi_ >& tk, const MagneticField* field, bool withErr) {
    Basic3DVector<float> pos(tk.POCA());
    GlobalPoint gpos(pos);
    Basic3DVector<float> mom(tk.momentum());
    GlobalVector gmom(mom);
    // no track charge so used track transverse curvature curvature
    int dummy=0;
    float speedOfLightConverted = CLHEP::c_light/1.0E5; // B*c/2E11 - converts q/pt to track angle at some radius from beamline
    double mMagneticFieldStrength = field->inTesla(GlobalPoint(0,0,0)).z();
    float trk_signedPt = speedOfLightConverted * mMagneticFieldStrength / tk.rInv(); // transverse curvature
    GlobalTrajectoryParameters par(gpos, gmom, trk_signedPt, dummy, field);
    //if (!withErr)
    //  return FreeTrajectoryState(par);
    //CurvilinearTrajectoryError newError;  // zeroed                                                                                                                            
    //auto& C = newError.matrix();
    AlgebraicSymMatrix55 mat;

    int dim = 5;
    int k = 0;
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j <= i; j++) {
	mat(i, j) = 1.0E-8;
      }
    }

    CurvilinearTrajectoryError newError(mat);

    // auto sin2Theta = par.momentum().perp2() / par.momentum().mag2(); 
    // AlgebraicSymMatrix55 mat;                                                                                                                                                  
    // double maxC00 = 1.0;
    // C[0][0] = std::max(0.1*sin2Theta, maxC00);
    // C[1][1] = std::abs(0.1*atan(tk.tanL()));
    // C[2][2] = std::abs(0.1*tk.phi());
    // C[3][3] = std::abs(0.1*tk.d0());
    // C[4][4] = std::abs(0.1*tk.z0());

    return FreeTrajectoryState(par, newError);
  }

  FreeTrajectoryState innerFreeState(const reco::Track& tk, const MagneticField* field, bool withErr) {
    Basic3DVector<float> pos(tk.innerPosition());
    GlobalPoint gpos(pos);
    Basic3DVector<float> mom(tk.innerMomentum());
    GlobalVector gmom(mom);
    GlobalTrajectoryParameters par(gpos, gmom, tk.charge(), field);
    if (!withErr)
      return FreeTrajectoryState(par);
    CurvilinearTrajectoryError err(tk.extra()->innerStateCovariance());
    return FreeTrajectoryState(par, err);
  }

  FreeTrajectoryState outerFreeState(const reco::Track& tk, const MagneticField* field, bool withErr) {
    Basic3DVector<float> pos(tk.outerPosition());
    GlobalPoint gpos(pos);
    Basic3DVector<float> mom(tk.outerMomentum());
    GlobalVector gmom(mom);
    GlobalTrajectoryParameters par(gpos, gmom, tk.charge(), field);
    if (!withErr)
      return FreeTrajectoryState(par);
    CurvilinearTrajectoryError err(tk.extra()->outerStateCovariance());
    return FreeTrajectoryState(par, err);
  }

  TrajectoryStateOnSurface innerStateOnSurface(const reco::Track& tk,
                                               const TrackingGeometry& geom,
                                               const MagneticField* field,
                                               bool withErr) {
    const Surface& surface = geom.idToDet(DetId(tk.extra()->innerDetId()))->surface();
    return TrajectoryStateOnSurface(innerFreeState(tk, field, withErr), surface);
  }

  TrajectoryStateOnSurface outerStateOnSurface(const reco::Track& tk,
                                               const TrackingGeometry& geom,
                                               const MagneticField* field,
                                               bool withErr) {
    const Surface& surface = geom.idToDet(DetId(tk.extra()->outerDetId()))->surface();
    return TrajectoryStateOnSurface(outerFreeState(tk, field, withErr), surface);
  }

}  // namespace trajectoryStateTransform
