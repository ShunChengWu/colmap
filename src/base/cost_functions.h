// Copyright (c) 2022, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#ifndef COLMAP_SRC_BASE_COST_FUNCTIONS_H_
#define COLMAP_SRC_BASE_COST_FUNCTIONS_H_

#include <Eigen/Core>

#include <ceres/ceres.h>
#include <ceres/rotation.h>


//#include <vector>
#include "base/intensityPatch.h"
//#include "optim/photometric_bundle_adjustment.h"  //TODO: maybe move this to a seperate file?
#include "util/pba.h"
#include "util/types.h"
//#include <Eigen/Geometry>

namespace colmap {

// Standard bundle adjustment cost function for variable
// camera pose and calibration and point parameters.
template <typename CameraModel>
class BundleAdjustmentCostFunction {
 public:
  explicit BundleAdjustmentCostFunction(const Eigen::Vector2d& point2D)
      : observed_x_(point2D(0)), observed_y_(point2D(1)) {}

  static ceres::CostFunction* Create(const Eigen::Vector2d& point2D) {
    return (new ceres::AutoDiffCostFunction<
            BundleAdjustmentCostFunction<CameraModel>, 2, 4, 3, 3,
            CameraModel::kNumParams>(
        new BundleAdjustmentCostFunction(point2D)));
  }

  template <typename T>
  bool operator()(const T* const qvec, const T* const tvec,
                  const T* const point3D, const T* const camera_params,
                  T* residuals) const {
    // Rotate and translate.
    T projection[3];
    ceres::UnitQuaternionRotatePoint(qvec, point3D, projection);
    projection[0] += tvec[0];
    projection[1] += tvec[1];
    projection[2] += tvec[2];

    // Project to image plane.
    projection[0] /= projection[2];
    projection[1] /= projection[2];

    // Distort and transform to pixel space.
    CameraModel::WorldToImage(camera_params, projection[0], projection[1],
                              &residuals[0], &residuals[1]);

    // Re-projection error.
    residuals[0] -= T(observed_x_);
    residuals[1] -= T(observed_y_);

    return true;
  }

 private:
  const double observed_x_;
  const double observed_y_;
};

// Bundle adjustment cost function for variable
// camera calibration and point parameters, and fixed camera pose.
template <typename CameraModel>
class BundleAdjustmentConstantPoseCostFunction {
 public:
  BundleAdjustmentConstantPoseCostFunction(const Eigen::Vector4d& qvec,
                                           const Eigen::Vector3d& tvec,
                                           const Eigen::Vector2d& point2D)
      : qw_(qvec(0)),
        qx_(qvec(1)),
        qy_(qvec(2)),
        qz_(qvec(3)),
        tx_(tvec(0)),
        ty_(tvec(1)),
        tz_(tvec(2)),
        observed_x_(point2D(0)),
        observed_y_(point2D(1)) {}

  static ceres::CostFunction* Create(const Eigen::Vector4d& qvec,
                                     const Eigen::Vector3d& tvec,
                                     const Eigen::Vector2d& point2D) {
    return (new ceres::AutoDiffCostFunction<
            BundleAdjustmentConstantPoseCostFunction<CameraModel>, 2, 3,
            CameraModel::kNumParams>(
        new BundleAdjustmentConstantPoseCostFunction(qvec, tvec, point2D)));
  }

  template <typename T>
  bool operator()(const T* const point3D, const T* const camera_params,
                  T* residuals) const {
    const T qvec[4] = {T(qw_), T(qx_), T(qy_), T(qz_)};

    // Rotate and translate.
    T projection[3];
    ceres::UnitQuaternionRotatePoint(qvec, point3D, projection);
    projection[0] += T(tx_);
    projection[1] += T(ty_);
    projection[2] += T(tz_);

    // Project to image plane.
    projection[0] /= projection[2];
    projection[1] /= projection[2];

    // Distort and transform to pixel space.
    CameraModel::WorldToImage(camera_params, projection[0], projection[1],
                              &residuals[0], &residuals[1]);

    // Re-projection error.
    residuals[0] -= T(observed_x_);
    residuals[1] -= T(observed_y_);

    return true;
  }

 private:
  const double qw_;
  const double qx_;
  const double qy_;
  const double qz_;
  const double tx_;
  const double ty_;
  const double tz_;
  const double observed_x_;
  const double observed_y_;
};


template<class CameraModel, class T> bool photometricLoss(
    const IntensityPatch &patch_i,
    const Eigen::Vector2i &boundary,
    const T* const qvec_i, const T* const tvec_i,
    const T* const qvec_j, const T* const tvec_j,
    const T* const camera_params_i, const T* const camera_params_j,
    T const* const ab_i, T const * const ab_j,
    const std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>> &interpolator_j,
    T const* const invD,
    bool useInvGradWeight, double weightC,
    T* sResiduals
){


  /*build affine*/
  Eigen::Matrix<T,2,1> affine;
  const T t_j = T(1), t_i = T(1);//assume exposure time are the same
  affine.x() = ceres::exp(ab_j[0]-ab_i[0])*t_j/t_i;
  affine.y() = ab_j[1] - affine.x()*ab_i[1];

  Eigen::Matrix<T,NumResidualPattern,1> I_j;
  T intensity;
  Eigen::Map<Eigen::Matrix<T, NumResidualPattern, 1>> residuals(sResiduals);
  Eigen::Map<Eigen::Quaternion<T> const> const eigen_qvec_i(qvec_i), eigen_qvec_j(qvec_j);
  Eigen::Map<Eigen::Vector3<T> const> const eigen_tvec_i(tvec_i), eigen_tvec_j(tvec_j);

  T depth = T(1)/ (*invD);
  for(size_t i=0;i<NumResidualPattern;++i) {
    /*reproject to world*/
    Eigen::Vector3<T> pt3D;
    pt3D.z() = T(1);

    CameraModel::ImageToWorld(
        camera_params_i,
        patch_i.coords.col(i).cast<T>()[0],
        patch_i.coords.col(i).cast<T>()[1],
        &pt3D[0],
        &pt3D[1]);
    pt3D *= depth;

    /*project from frame coord. i to j*/
//    // from i to world
    pt3D = eigen_qvec_i.inverse() * pt3D;
    pt3D += -(eigen_qvec_i.inverse()*eigen_tvec_i);
    // from world to i
    pt3D = eigen_qvec_j * pt3D;
    pt3D += eigen_tvec_j;

//    pt3D = T_j_i.template topLeftCorner<3,3>() * pt3D + T_j_i.template topRightCorner<3,1>();

    /*project to image*/
    pt3D[0] /= pt3D[2];
    pt3D[1] /= pt3D[2];

    Eigen::Matrix<T,2,1> p2d_i_in_j;
    CameraModel::WorldToImage(camera_params_j, pt3D[0],pt3D[1],
                              &p2d_i_in_j[0],&p2d_i_in_j[1]);

    interpolator_j->Evaluate(p2d_i_in_j.y(), p2d_i_in_j.x(), &intensity);
    I_j.row(i) << intensity;
  }

  /*calculate reverse weight according to DSO*/
  Eigen::Matrix<T,NumResidualPattern,1> Weight = Eigen::Matrix<T,NumResidualPattern,1>::Ones();
  if(useInvGradWeight)
  {
    T C2 = T(weightC);
    auto absGrads = patch_i.gradients.cwiseAbs().cast<T>() / T(255);
    Weight.array() = C2 / (C2+absGrads.array()*absGrads.array());
  }

  residuals = Weight.array() * (I_j.array()-(affine.x()*patch_i.intensities.array()+affine.y()));
//  for (size_t i=0;i<NumResidualPattern;++i) {
//    residuals.row(i) = Weight.row(i) * residuals[i];
//  }
  return true;
}

// Photometric bundle adjustment cost function
template <typename CameraModel>
class PhotometricBundleAdjustmentCostFunction {
 public:
  explicit PhotometricBundleAdjustmentCostFunction(
      const IntensityPatch & patch_i,
      std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>> interpolator_j,
      Eigen::Vector2i imageBoundary, double constC)
      : patch_i_(patch_i),
        imageBoundary_(std::move(imageBoundary)),
        interpolator_j_(std::move(interpolator_j)),
        constC_(constC)
  {}

  static ceres::CostFunction* CreateOneCam(const IntensityPatch &patch_i,
                                     std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>> interpolator_j,
                                     Eigen::Vector2i imageBoundary,
                                           double constC) {
    return (new ceres::AutoDiffCostFunction<
            PhotometricBundleAdjustmentCostFunction<CameraModel>, NumResidualPattern,
                4, 3, 4, 3, CameraModel::kNumParams,
            2,2,1>(
        new PhotometricBundleAdjustmentCostFunction(
            patch_i,interpolator_j,imageBoundary,constC)));
  }

  static ceres::CostFunction* CreateTwoCams(const IntensityPatch &patch_i,
                                     std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>> interpolator_j,
                                     Eigen::Vector2i imageBoundary,
                                            double constC) {
    return (new ceres::AutoDiffCostFunction<
            PhotometricBundleAdjustmentCostFunction<CameraModel>, NumResidualPattern,
            4, 3, 4, 3, CameraModel::kNumParams,CameraModel::kNumParams,
            2,2,1>(
        new PhotometricBundleAdjustmentCostFunction(
            patch_i,interpolator_j,imageBoundary,constC)));
  }

  /*One camera*/
  template <typename T>
  bool operator()(const T* const qvec_i, const T* const tvec_i,
                  const T* const qvec_j, const T* const tvec_j,
                  const T* const camera_params,
                  const T* const ab_i, const T* const ab_j,
                  const T* const invserseD,
                  T* residuals) const {
    return photometricLoss<CameraModel,T>(
        patch_i_,imageBoundary_,
        qvec_i,tvec_i,
        qvec_j,tvec_j,
        camera_params,camera_params,
        ab_i,ab_j,interpolator_j_,invserseD,true,constC_,residuals
    );
  }

  /*Two Cam*/
  template <typename T>
  bool operator()(const T* const qvec_i, const T* const tvec_i,
                  const T* const qvec_j, const T* const tvec_j,
                  const T* const camera_params_i, const T* const camera_params_j,
                  const T* const ab_i, const T* const ab_j,
                  const T* const invserseD,
                  T* residuals) const {
    return photometricLoss<CameraModel,T>(
        patch_i_,imageBoundary_,
        qvec_i,tvec_i,
        qvec_j,tvec_j,
        camera_params_i,camera_params_j,
        ab_i,ab_j,interpolator_j_,invserseD,true,constC_,residuals
    );
  }

 private:
  IntensityPatch patch_i_;
  Eigen::Vector2i imageBoundary_;
  std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>>
      interpolator_j_;
  double constC_;
};

template <typename CameraModel>
class PhotometricBundleAdjustmentConstantSrcCostFunction {
 public:
  explicit PhotometricBundleAdjustmentConstantSrcCostFunction(
      const Eigen::Vector4d& qvec_i,
      const Eigen::Vector3d& tvec_i,
      const IntensityPatch & patch_i,
      std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>> interpolator_j,
      Eigen::Vector2i imageBoundary,
      double constC)
      : qw_i_(qvec_i(0)),
        qx_i_(qvec_i(1)),
        qy_i_(qvec_i(2)),
        qz_i_(qvec_i(3)),
        tx_i_(tvec_i(0)),
        ty_i_(tvec_i(1)),
        tz_i_(tvec_i(2)),
        patch_i_(patch_i),
        imageBoundary_(std::move(imageBoundary)),
        interpolator_j_(std::move(interpolator_j)),
        constC_(constC)
  {}

  static ceres::CostFunction* CreateOneCam(const Eigen::Vector4d& qvec_i,
                                     const Eigen::Vector3d& tvec_i,
                                     const IntensityPatch &patch_i,
                                     std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>> interpolator_j,
                                     Eigen::Vector2i imageBoundary,
                                           double constC) {
    return (new ceres::AutoDiffCostFunction<
            PhotometricBundleAdjustmentConstantSrcCostFunction<CameraModel>, NumResidualPattern,
            4, 3, CameraModel::kNumParams,
            2,2,1>(
        new PhotometricBundleAdjustmentConstantSrcCostFunction(
            qvec_i,tvec_i,patch_i,interpolator_j,imageBoundary,constC)));
  }

  static ceres::CostFunction* CreateTwoCams(const Eigen::Vector4d& qvec_i,
                                     const Eigen::Vector3d& tvec_i,
                                     const IntensityPatch &patch_i,
                                     std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>> interpolator_j,
                                     Eigen::Vector2i imageBoundary,
                                            double constC) {
    return (new ceres::AutoDiffCostFunction<
            PhotometricBundleAdjustmentConstantSrcCostFunction<CameraModel>, NumResidualPattern,
            4, 3, CameraModel::kNumParams,CameraModel::kNumParams,
            2,2,1>(
        new PhotometricBundleAdjustmentConstantSrcCostFunction(
            qvec_i,tvec_i,patch_i,interpolator_j,imageBoundary,constC)));
  }

  /*One camera*/
  template <typename T>
  bool operator()(const T* const qvec_j, const T* const tvec_j,
                  const T* const camera_params,
                  const T* const ab_i, const T* const ab_j,
                  const T* const invserseD,
                  T* residuals) const {
    const T qvec_i[4] = {T(qw_i_), T(qx_i_), T(qy_i_), T(qz_i_)};
    const T tvec_i[3] = {T(tx_i_), T(ty_i_), T(tz_i_)};

    return photometricLoss<CameraModel,T>(
        patch_i_,imageBoundary_,
        qvec_i,tvec_i,
        qvec_j,tvec_j,
        camera_params,camera_params,
        ab_i,ab_j,interpolator_j_,invserseD,true,constC_,residuals
    );
  }

  /*Two cameras*/
  template <typename T>
  bool operator()(const T* const qvec_j, const T* const tvec_j,
                  const T* const camera_params_i, const T* const camera_params_j,
                  const T* const ab_i, const T* const ab_j,
                  const T* const invserseD,
                  T* residuals) const {
    const T qvec_i[4] = {T(qw_i_), T(qx_i_), T(qy_i_), T(qz_i_)};
    const T tvec_i[3] = {T(tx_i_), T(ty_i_), T(tz_i_)};

    return photometricLoss<CameraModel,T>(
        patch_i_,imageBoundary_,
        qvec_i,tvec_i,
        qvec_j,tvec_j,
        camera_params_i,camera_params_j,
        ab_i,ab_j,interpolator_j_,invserseD,true,constC_,residuals
    );
  }

 private:
  const double qw_i_;
  const double qx_i_;
  const double qy_i_;
  const double qz_i_;
  const double tx_i_;
  const double ty_i_;
  const double tz_i_;
  IntensityPatch patch_i_;
  Eigen::Vector2i imageBoundary_;
  std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>>
      interpolator_j_;
  double constC_;
};

template <typename CameraModel>
class PhotometricBundleAdjustmentConstantTgtCostFunction {
 public:
  explicit PhotometricBundleAdjustmentConstantTgtCostFunction(
      const Eigen::Vector4d& qvec_j,
      const Eigen::Vector3d& tvec_j,
      const IntensityPatch & patch_i,
      std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>> interpolator_j,
      Eigen::Vector2i imageBoundary,
      double constC)
      : qw_j_(qvec_j(0)),
        qx_j_(qvec_j(1)),
        qy_j_(qvec_j(2)),
        qz_j_(qvec_j(3)),
        tx_j_(tvec_j(0)),
        ty_j_(tvec_j(1)),
        tz_j_(tvec_j(2)),
        patch_i_(patch_i),
        imageBoundary_(std::move(imageBoundary)),
        interpolator_j_(std::move(interpolator_j)),
        constC_(constC)
  {}

  static ceres::CostFunction* CreateOneCam(const Eigen::Vector4d& qvec_j,
                                     const Eigen::Vector3d& tvec_j,
                                     const IntensityPatch &patch_i,
                                     std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>> interpolator_j,
                                     Eigen::Vector2i imageBoundary,
                                           double constC) {
    return (new ceres::AutoDiffCostFunction<
            PhotometricBundleAdjustmentConstantTgtCostFunction<CameraModel>, NumResidualPattern,
            4, 3, CameraModel::kNumParams,
            2,2,1>(
        new PhotometricBundleAdjustmentConstantTgtCostFunction(
            qvec_j,tvec_j,patch_i,interpolator_j,imageBoundary,constC)));
  }

  static ceres::CostFunction* CreateTwoCams(const Eigen::Vector4d& qvec_j,
                                     const Eigen::Vector3d& tvec_j,
                                     const IntensityPatch &patch_i,
                                     std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>> interpolator_j,
                                     Eigen::Vector2i imageBoundary,
                                            double constC) {
    return (new ceres::AutoDiffCostFunction<
            PhotometricBundleAdjustmentConstantTgtCostFunction<CameraModel>, NumResidualPattern,
            4, 3, CameraModel::kNumParams,CameraModel::kNumParams,
            2,2,1>(
        new PhotometricBundleAdjustmentConstantTgtCostFunction(
            qvec_j,tvec_j,patch_i,interpolator_j,imageBoundary,constC)));
  }

  /*One camera*/
  template <typename T>
  bool operator()(const T* const qvec_i, const T* const tvec_i,
                  const T* const camera_params,
                  const T* const ab_i, const T* const ab_j,
                  const T* const invserseD,
                  T* residuals) const {
    const T qvec_j[4] = {T(qw_j_), T(qx_j_), T(qy_j_), T(qz_j_)};
    const T tvec_j[3] = {T(tx_j_), T(ty_j_), T(tz_j_)};

    return photometricLoss<CameraModel,T>(
        patch_i_,imageBoundary_,
        qvec_i,tvec_i,
        qvec_j,tvec_j,
        camera_params,camera_params,
        ab_i,ab_j,interpolator_j_,invserseD,true,constC_,residuals
    );
  }

  /*Two Cams*/
  template <typename T>
  bool operator()(const T* const qvec_i, const T* const tvec_i,
                  const T* const camera_params_i, const T* const camera_params_j,
                  const T* const ab_i, const T* const ab_j,
                  const T* const invserseD,
                  T* residuals) const {
    const T qvec_j[4] = {T(qw_j_), T(qx_j_), T(qy_j_), T(qz_j_)};
    const T tvec_j[3] = {T(tx_j_), T(ty_j_), T(tz_j_)};

    return photometricLoss<CameraModel,T>(
        patch_i_,imageBoundary_,
        qvec_i,tvec_i,
        qvec_j,tvec_j,
        camera_params_i,camera_params_j,
        ab_i,ab_j,interpolator_j_,invserseD,true,constC_,residuals
    );
  }

 private:
  const double qw_j_;
  const double qx_j_;
  const double qy_j_;
  const double qz_j_;
  const double tx_j_;
  const double ty_j_;
  const double tz_j_;
  IntensityPatch patch_i_;
  Eigen::Vector2i imageBoundary_;
  std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>>
      interpolator_j_;
  double constC_;
};

template <typename CameraModel>
class PhotometricBundleAdjustmentConstantCostFunction {
 public:
  explicit PhotometricBundleAdjustmentConstantCostFunction(
      const Eigen::Vector4d& qvec_i,
      const Eigen::Vector3d& tvec_i,
      const Eigen::Vector4d& qvec_j,
      const Eigen::Vector3d& tvec_j,
      const IntensityPatch & patch_i,
      std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>> interpolator_j,
      Eigen::Vector2i imageBoundary,
      double constC)
      : qw_i_(qvec_i(0)),
        qx_i_(qvec_i(1)),
        qy_i_(qvec_i(2)),
        qz_i_(qvec_i(3)),
        tx_i_(tvec_i(0)),
        ty_i_(tvec_i(1)),
        tz_i_(tvec_i(2)),
        qw_j_(qvec_j(0)),
        qx_j_(qvec_j(1)),
        qy_j_(qvec_j(2)),
        qz_j_(qvec_j(3)),
        tx_j_(tvec_j(0)),
        ty_j_(tvec_j(1)),
        tz_j_(tvec_j(2)),
        patch_i_(patch_i),
        imageBoundary_(std::move(imageBoundary)),
        interpolator_j_(std::move(interpolator_j)),
        constC_(constC)
  {}

  static ceres::CostFunction* CreateOneCam(const Eigen::Vector4d& qvec_i,
                                     const Eigen::Vector3d& tvec_i,
                                     const Eigen::Vector4d& qvec_j,
                                     const Eigen::Vector3d& tvec_j,
                                     const IntensityPatch &patch_i,
                                     std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>> interpolator_j,
                                     Eigen::Vector2i imageBoundary,
                                           double constC) {
    return (new ceres::AutoDiffCostFunction<
            PhotometricBundleAdjustmentConstantCostFunction<CameraModel>, NumResidualPattern,
            CameraModel::kNumParams,
            2,2,1>(
        new PhotometricBundleAdjustmentConstantCostFunction(
            qvec_i,tvec_i,qvec_j,tvec_j,patch_i,interpolator_j,imageBoundary,constC)));
  }

  static ceres::CostFunction* CreateTwoCams(const Eigen::Vector4d& qvec_i,
                                     const Eigen::Vector3d& tvec_i,
                                     const Eigen::Vector4d& qvec_j,
                                     const Eigen::Vector3d& tvec_j,
                                     const IntensityPatch &patch_i,
                                     std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>> interpolator_j,
                                     Eigen::Vector2i imageBoundary,
                                            double constC) {
    return (new ceres::AutoDiffCostFunction<
            PhotometricBundleAdjustmentConstantCostFunction<CameraModel>, NumResidualPattern,
            CameraModel::kNumParams,CameraModel::kNumParams,
            2,2,1>(
        new PhotometricBundleAdjustmentConstantCostFunction(
            qvec_i,tvec_i,qvec_j,tvec_j,patch_i,interpolator_j,imageBoundary,constC)));
  }

  /*One camera*/
  template <typename T>
  bool operator()(
                  const T* const camera_params,
                  const T* const ab_i, const T* const ab_j,
                  const T* const invserseD,
                  T* residuals) const {
    const T qvec_i[4] = {T(qw_i_), T(qx_i_), T(qy_i_), T(qz_i_)};
    const T tvec_i[3] = {T(tx_i_), T(ty_i_), T(tz_i_)};
    const T qvec_j[4] = {T(qw_j_), T(qx_j_), T(qy_j_), T(qz_j_)};
    const T tvec_j[3] = {T(tx_j_), T(ty_j_), T(tz_j_)};

    return photometricLoss<CameraModel,T>(
        patch_i_,imageBoundary_,
        qvec_i,tvec_i,
        qvec_j,tvec_j,
        camera_params,camera_params,
        ab_i,ab_j,interpolator_j_,invserseD,true,constC_,residuals
    );
  }

  /*Two cameras*/
  template <typename T>
  bool operator()(
      const T* const camera_params_i, const T* const camera_params_j,
      const T* const ab_i, const T* const ab_j,
      const T* const invserseD,
      T* residuals) const {
    const T qvec_i[4] = {T(qw_i_), T(qx_i_), T(qy_i_), T(qz_i_)};
    const T tvec_i[3] = {T(tx_i_), T(ty_i_), T(tz_i_)};
    const T qvec_j[4] = {T(qw_j_), T(qx_j_), T(qy_j_), T(qz_j_)};
    const T tvec_j[3] = {T(tx_j_), T(ty_j_), T(tz_j_)};

    return photometricLoss<CameraModel,T>(
        patch_i_,imageBoundary_,
        qvec_i,tvec_i,
        qvec_j,tvec_j,
        camera_params_i,camera_params_j,
        ab_i,ab_j,interpolator_j_,invserseD,true,constC_,residuals
    );
  }

 private:
  const double qw_i_,qw_j_;
  const double qx_i_,qx_j_;
  const double qy_i_,qy_j_;
  const double qz_i_,qz_j_;
  const double tx_i_,tx_j_;
  const double ty_i_,ty_j_;
  const double tz_i_,tz_j_;
  IntensityPatch patch_i_;
  Eigen::Vector2i imageBoundary_;
  std::shared_ptr<ceres::BiCubicInterpolator<ceres::Grid2D<IntensityType, 1>>>
      interpolator_j_;
  double constC_;
};

// Rig bundle adjustment cost function for variable camera pose and calibration
// and point parameters. Different from the standard bundle adjustment function,
// this cost function is suitable for camera rigs with consistent relative poses
// of the cameras within the rig. The cost function first projects points into
// the local system of the camera rig and then into the local system of the
// camera within the rig.
template <typename CameraModel>
class RigBundleAdjustmentCostFunction {
 public:
  explicit RigBundleAdjustmentCostFunction(const Eigen::Vector2d& point2D)
      : observed_x_(point2D(0)), observed_y_(point2D(1)) {}

  static ceres::CostFunction* Create(const Eigen::Vector2d& point2D) {
    return (new ceres::AutoDiffCostFunction<
            RigBundleAdjustmentCostFunction<CameraModel>, 2, 4, 3, 4, 3, 3,
            CameraModel::kNumParams>(
        new RigBundleAdjustmentCostFunction(point2D)));
  }

  template <typename T>
  bool operator()(const T* const rig_qvec, const T* const rig_tvec,
                  const T* const rel_qvec, const T* const rel_tvec,
                  const T* const point3D, const T* const camera_params,
                  T* residuals) const {
    // Concatenate rotations.
    T qvec[4];
    ceres::QuaternionProduct(rel_qvec, rig_qvec, qvec);

    // Concatenate translations.
    T tvec[3];
    ceres::UnitQuaternionRotatePoint(rel_qvec, rig_tvec, tvec);
    tvec[0] += rel_tvec[0];
    tvec[1] += rel_tvec[1];
    tvec[2] += rel_tvec[2];

    // Rotate and translate.
    T projection[3];
    ceres::UnitQuaternionRotatePoint(qvec, point3D, projection);
    projection[0] += tvec[0];
    projection[1] += tvec[1];
    projection[2] += tvec[2];

    // Project to image plane.
    projection[0] /= projection[2];
    projection[1] /= projection[2];

    // Distort and transform to pixel space.
    CameraModel::WorldToImage(camera_params, projection[0], projection[1],
                              &residuals[0], &residuals[1]);

    // Re-projection error.
    residuals[0] -= T(observed_x_);
    residuals[1] -= T(observed_y_);

    return true;
  }

 private:
  const double observed_x_;
  const double observed_y_;
};

// Cost function for refining two-view geometry based on the Sampson-Error.
//
// First pose is assumed to be located at the origin with 0 rotation. Second
// pose is assumed to be on the unit sphere around the first pose, i.e. the
// pose of the second camera is parameterized by a 3D rotation and a
// 3D translation with unit norm. `tvec` is therefore over-parameterized as is
// and should be down-projected using `SphereManifold`.
class RelativePoseCostFunction {
 public:
  RelativePoseCostFunction(const Eigen::Vector2d& x1, const Eigen::Vector2d& x2)
      : x1_(x1(0)), y1_(x1(1)), x2_(x2(0)), y2_(x2(1)) {}

  static ceres::CostFunction* Create(const Eigen::Vector2d& x1,
                                     const Eigen::Vector2d& x2) {
    return (new ceres::AutoDiffCostFunction<RelativePoseCostFunction, 1, 4, 3>(
        new RelativePoseCostFunction(x1, x2)));
  }

  template <typename T>
  bool operator()(const T* const qvec, const T* const tvec,
                  T* residuals) const {
    Eigen::Matrix<T, 3, 3, Eigen::RowMajor> R;
    ceres::QuaternionToRotation(qvec, R.data());

    // Matrix representation of the cross product t x R.
    Eigen::Matrix<T, 3, 3> t_x;
    t_x << T(0), -tvec[2], tvec[1], tvec[2], T(0), -tvec[0], -tvec[1], tvec[0],
        T(0);

    // Essential matrix.
    const Eigen::Matrix<T, 3, 3> E = t_x * R;

    // Homogeneous image coordinates.
    const Eigen::Matrix<T, 3, 1> x1_h(T(x1_), T(y1_), T(1));
    const Eigen::Matrix<T, 3, 1> x2_h(T(x2_), T(y2_), T(1));

    // Squared sampson error.
    const Eigen::Matrix<T, 3, 1> Ex1 = E * x1_h;
    const Eigen::Matrix<T, 3, 1> Etx2 = E.transpose() * x2_h;
    const T x2tEx1 = x2_h.transpose() * Ex1;
    residuals[0] = x2tEx1 * x2tEx1 /
                   (Ex1(0) * Ex1(0) + Ex1(1) * Ex1(1) + Etx2(0) * Etx2(0) +
                    Etx2(1) * Etx2(1));

    return true;
  }

 private:
  const double x1_;
  const double y1_;
  const double x2_;
  const double y2_;
};

inline void SetQuaternionManifold(ceres::Problem* problem, double* qvec) {
#if CERES_VERSION_MAJOR >= 2 && CERES_VERSION_MINOR >= 1
  problem->SetManifold(qvec, new ceres::QuaternionManifold);
#else
  problem->SetParameterization(qvec, new ceres::QuaternionParameterization);
#endif
}

inline void SetSubsetManifold(int size, const std::vector<int>& constant_params,
                              ceres::Problem* problem, double* params) {
#if CERES_VERSION_MAJOR >= 2 && CERES_VERSION_MINOR >= 1
  problem->SetManifold(params,
                       new ceres::SubsetManifold(size, constant_params));
#else
  problem->SetParameterization(
      params, new ceres::SubsetParameterization(size, constant_params));
#endif
}

template <int size>
inline void SetSphereManifold(ceres::Problem* problem, double* params) {
#if CERES_VERSION_MAJOR >= 2 && CERES_VERSION_MINOR >= 1
  problem->SetManifold(params, new ceres::SphereManifold<size>);
#else
  problem->SetParameterization(
      params, new ceres::HomogeneousVectorParameterization(size));
#endif
}

}  // namespace colmap

#endif  // COLMAP_SRC_BASE_COST_FUNCTIONS_H_
