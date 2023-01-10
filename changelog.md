# TODO
- [x] split image set when mapping
  - `--image_list_path` will do
- [x] Accepting known poses in estimation
- [ ] Combineing known and unknown pose
  - This should do with incremental mapping
- [x] Allow fixing existing image
  - This has provided API `fix_existing_images` in `IncrementalMapperOptions`


# Initial Pose Estimation
## Estimate initial pose pair
`controllers/incremental_mapper.cc:l431: const bool find_init_success = mapper.FindInitialImagePair(
            init_mapper_options, &image_id1, &image_id2);`
## Poses are estimated here
`sfm/incremental_mapper.cc:l189: EstimateInitialTwoViewGeometry`
## The estimated pose will be stored at 
`sfm/incremental_mapper.cc:l1195: prev_init_two_view_geometry_ = two_view_geometry;`

## Then the pose will be assigned to images here
`controllers/incremental_mapper.cc:l454: const bool reg_init_success = mapper.RegisterInitialImagePair`
`sfm/incremental_mapper.cc:l291: image2.Qvec() = prev_init_two_view_geometry_.qvec;`

Steps:
1. Allow IncrementalMapper estimate initial pair using given poses.
2. Allow Incremental mapping takes known images with known poses 
Try to make `sfm/incremental_mapper.cc:l1163:EstimateInitialTwoViewGeometry` works with known given pose
3. How to allow `IncrementalMapper::RegisterNextImage` taking known pose?
   - In `srm/incremental_mapper.cc:l519`, Continue tracks needs inliers to be added to `reconstruction_` and `triangulator_`.

# Logs
2023/1/9: Can do incremental mapping + triangulation without BA refinement. If BA is used, the scale will be wrong. 
I guess the problem is the no poses are fixed. Need to fix at least two to ensure correct scale

todo: check if point_triangulation fixes all poses, or does not optimize pose.

2023/1/10: Add fix given pose and all extrinsics option in mapper.
