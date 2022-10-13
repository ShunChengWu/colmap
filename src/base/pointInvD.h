//
// Created by sc on 10/9/22.
//

#ifndef COLMAP_SRC_BASE_INVERSEDEPTHPOINT_H_
#define COLMAP_SRC_BASE_INVERSEDEPTHPOINT_H_

#include <utility>
#include <vector>

#include <Eigen/Core>

#include "base/track.h"
#include "base/intensityPatch.h"
#include "util/logging.h"
#include "util/types.h"

namespace colmap {
  class PointInvD {
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    PointInvD();

    // The point coordinate in world space.
//    inline const Eigen::Vector3d& XYZ() const;
//    inline Eigen::Vector3d& XYZ();
//    inline double XYZ(const size_t idx) const;
//    inline double& XYZ(const size_t idx);
//    inline double X() const;
//    inline double Y() const;
//    inline double Z() const;
//    inline void SetXYZ(const Eigen::Vector3d& xyz);
    inline void SetInverseDepth(double invD);
    inline void SetHostFrame(TrackElement trackElement);
    inline void SetIntensityPatch(IntensityPatch intensityPatch);
    inline void SetTrack(const Track &track);

    inline double& InverseDepth();
    inline const double& InverseDepth() const;
//    inline double* InverseDepthPtr();

    inline TrackElement& HostFrame();
    inline const TrackElement& HostFrame() const;

    inline class IntensityPatch &IntensityPatch();
    inline const class IntensityPatch &IntensityPatch() const;

    inline class Track &Track();
    inline const class Track &Track() const;
//
//    // The RGB color of the point.
//    inline const Eigen::Vector3ub& Color() const;
//    inline Eigen::Vector3ub& Color();
//    inline uint8_t Color(const size_t idx) const;
//    inline uint8_t& Color(const size_t idx);
//    inline void SetColor(const Eigen::Vector3ub& color);
//
//    // The mean reprojection error in image space.
//    inline double Error() const;
//    inline bool HasError() const;
//    inline void SetError(const double error);
//
//    inline const class Track& Track() const;
//    inline class Track& Track();
//    inline void SetTrack(const class Track& track);

   private:
    // Inverse depth
    double inverseDepth_;

    // Whether this point is created from a geometric landmark
    bool isFromGeo_;

    // Host frame id
    TrackElement hostFrame_;

    // The pixel patches from the host frame
    class IntensityPatch intensityPatch_;

    // The track of the point has been observed by a list of images.
    class Track track_;
  };

  void PointInvD::SetInverseDepth(double invD) {inverseDepth_ = invD; }

  void PointInvD::SetHostFrame(TrackElement trackElement) {hostFrame_ = trackElement;}

  void PointInvD::SetTrack(const class Track &track) {track_ = track;}

  void PointInvD::SetIntensityPatch(class IntensityPatch intensityPatch) {intensityPatch_ = std::move(intensityPatch);};

  double& PointInvD::InverseDepth() {return inverseDepth_;}

//  double* PointInvD::InverseDepthPtr() {return &inverseDepth_;}

  const double& PointInvD::InverseDepth() const {return inverseDepth_;}

  TrackElement& PointInvD::HostFrame() {return hostFrame_;}

  const TrackElement& PointInvD::HostFrame() const {return hostFrame_;}

  class IntensityPatch &PointInvD::IntensityPatch() {return intensityPatch_;}

  const class IntensityPatch &PointInvD::IntensityPatch() const {return intensityPatch_;}

  class Track &PointInvD::Track() {return track_;}

  const Track &PointInvD::Track() const {return track_;}
}

#endif  // COLMAP_SRC_BASE_INVERSEDEPTHPOINT_H_
