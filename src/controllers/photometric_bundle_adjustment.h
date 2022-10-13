//
// Created by sc on 10/9/22.
//

#ifndef COLMAP_SRC_CONTROLLERS_PHOTOMETRIC_BUNDLE_ADJUSTMENT_H
#define COLMAP_SRC_CONTROLLERS_PHOTOMETRIC_BUNDLE_ADJUSTMENT_H

#include "base/reconstruction.h"
#include "util/option_manager.h"
#include "util/threading.h"

namespace colmap {

// Class that controls the global bundle adjustment procedure.
class PhotometricBundleAdjustmentController : public Thread {
 public:
  PhotometricBundleAdjustmentController(const OptionManager& options,
                             Reconstruction* reconstruction);

 private:
  void Run();

  const OptionManager options_;
  Reconstruction* reconstruction_;
};

}  // namespace colmap

#endif  // COLMAP_SRC_CONTROLLERS_PHOTOMETRIC_BUNDLE_ADJUSTMENT_H
