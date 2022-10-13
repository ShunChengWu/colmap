//
// Created by sc on 10/9/22.
//

#ifndef COLMAP_SRC_BASE_INTENSITYPATCH_H_
#define COLMAP_SRC_BASE_INTENSITYPATCH_H_

#include "util/types.h"
#include <Eigen/Core>

namespace colmap {

  struct IntensityPatch {
    PatternIntensities intensities;
    Gradients gradients;
    PatternCoords coords;
  };

}

#endif  // COLMAP_SRC_BASE_INTENSITYPATCH_H_
