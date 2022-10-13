//
// Created by sc on 10/10/22.
//

#ifndef COLMAP_UTIL_PBA_H
#define COLMAP_UTIL_PBA_H
#include <Eigen/Core>
#include "util/types.h"
#include "base/image.h"
#include "base/intensityPatch.h"
namespace colmap {
[[maybe_unused]] static const Eigen::Vector2i
    PhotometricResidualPatterns[NumResidualPattern] = {
        Eigen::Vector2i(0, 0),  Eigen::Vector2i(2, 0),   Eigen::Vector2i(1, -1),
        Eigen::Vector2i(0, -2), Eigen::Vector2i(-1, -1), Eigen::Vector2i(-2, 0),
        Eigen::Vector2i(-1, 1), Eigen::Vector2i(0, 2),
};
static bool GetPatchInfo(const Image& image, const Point2D& p2d,
                         IntensityPatch& patch) {
  IntensityType intensity, gradient;


  for (size_t i = 0; i < NumResidualPattern; ++i) {
    auto p2d_i = p2d.XY() + PhotometricResidualPatterns[i].cast<double>();

    // get intensity
    image.IntensityInterpolator()->Evaluate(p2d_i.y(),p2d_i.x(), &intensity);

    // get gradient
    image.IntensityGradientInterpolator()->Evaluate(p2d_i.y(),p2d_i.x(),&gradient);

    patch.intensities.row(i) << intensity;
    patch.coords.col(i) = p2d_i;
    patch.gradients.row(i) << gradient;
  }

  // evaluate the mean gradient, skip if it's too large
  // the intuition here is similar to the inverse grad weighting in the DSO paper
  // if the gradient is too larger a point may locate at the boundary, which
  // may be unstable for PBA
  auto meanGrad = patch.gradients.cwiseAbs().sum()/IntensityType (255);
  if(meanGrad > 0.9 ) return false;
  return true;
}
}

#endif  // COLMAP_UTIL_PBA_H
