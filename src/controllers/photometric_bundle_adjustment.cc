//
// Created by sc on 10/9/22.
//

#include "controllers/photometric_bundle_adjustment.h"

#include <ceres/ceres.h>

#include "optim/photometric_bundle_adjustment.h"
#include "util/misc.h"

using namespace colmap;// {

namespace {

// Callback functor called after each bundle adjustment iteration.
class PhotometricBundleAdjustmentIterationCallback
    : public ceres::IterationCallback {
 public:
  explicit PhotometricBundleAdjustmentIterationCallback(Thread* thread)
      : thread_(thread) {}

  virtual ceres::CallbackReturnType operator()(
      const ceres::IterationSummary& summary) {
    CHECK_NOTNULL(thread_);
    thread_->BlockIfPaused();
    if (thread_->IsStopped()) {
      return ceres::SOLVER_TERMINATE_SUCCESSFULLY;
    } else {
      return ceres::SOLVER_CONTINUE;
    }
  }

 private:
  Thread* thread_;
};
}

PhotometricBundleAdjustmentController::PhotometricBundleAdjustmentController(
    const OptionManager& options, Reconstruction* reconstruction)
    : options_(options), reconstruction_(reconstruction) {}

void PhotometricBundleAdjustmentController::Run() {
  CHECK_NOTNULL(reconstruction_);

  PrintHeading1("Global photometric bundle adjustment");

  const std::vector<image_t>& reg_image_ids = reconstruction_->RegImageIds();

  if (reg_image_ids.size() < 2) {
    std::cout << "ERROR: Need at least two views." << std::endl;
    return;
  }

  // Avoid degeneracies in bundle adjustment.
  reconstruction_->FilterObservationsWithNegativeDepth();

  PhotometricBundleAdjustmentOptions pba_options = *options_.photometric_bundle_adjustment;
  pba_options.solver_options.minimizer_progress_to_stdout = true;

  PhotometricBundleAdjustmentIterationCallback iteration_callback(this);
  pba_options.solver_options.callbacks.push_back(&iteration_callback);

  // Configure photometric bundle adjustment.
  PhotometricBundleAdjustmentConfig pba_config;
  for (const image_t image_id : reg_image_ids) {
    pba_config.AddImage(image_id);
  }

  if (pba_options.fixed_poses.empty()) {
    pba_config.SetConstantPose(reg_image_ids[0]);
    pba_config.SetConstantTvec(reg_image_ids[1], {0});
  } else {
    std::stringstream ss;
    for (auto idx : pba_options.fixed_poses) ss << idx << " ";
    LOG(INFO) << "use custom fixed poses with image ids: " << ss.str();
    for (auto idx : pba_options.fixed_poses) pba_config.SetConstantPose(idx);
  }

  // Run bundle adjustment.
  PhotometricBundleAdjuster photometric_bundle_adjuster(pba_options,
                                                        pba_config);
  photometric_bundle_adjuster.Solve(reconstruction_);

  GetTimer().PrintMinutes();
}

//}