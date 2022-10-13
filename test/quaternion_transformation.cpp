//
// Created by sc on 10/13/22.
//
#include <Eigen/Geometry>
#include <random>
#include <iostream>
#include <iostream>

int main(){
  // generate random rotation and translation
  std::random_device dev;
  std::mt19937 gen(dev());
  std::uniform_real_distribution<> dist(0.0,1.0); // distribution in range [1, 6]

  Eigen::Quaterniond quat(dist(gen),dist(gen),dist(gen),dist(gen));
  quat.normalize();
  Eigen::Vector3d translation(dist(gen),dist(gen),dist(gen));

  /**/
  Eigen::Vector3d pt (dist(gen),dist(gen),dist(gen));

  // forward transform
  Eigen::Vector3d pt_f = quat*pt + translation;

  // backward transform
  Eigen::Vector3d pt_b_r = quat.inverse()*pt_f;
  Eigen::Vector3d pt_b_t = -(quat.inverse()*translation);
  Eigen::Vector3d pt_b = pt_b_r + pt_b_t;

  // residual
  double diff = (pt-pt_b).norm();


  /*compare with SE3*/
  Eigen::Matrix4d mat44 = Eigen::Matrix4d::Identity();
  mat44.topLeftCorner<3,3>() = quat.toRotationMatrix();
  mat44.topRightCorner<3,1>() = translation;

  // forward
  Eigen::Vector3d pt_f_ = (mat44 * pt.homogeneous()).head<3>();
//  Eigen::Vector3d pt_f_ = mat44.topLeftCorner<3,3>() * pt + mat44.topRightCorner<3,1>();

  // backward
  Eigen::Matrix4d mat44_inv = mat44.inverse();
  Eigen::Vector3d pt_b_r_ = mat44_inv.topLeftCorner<3,3>() * pt_f_;
  Eigen::Vector3d pt_b_t_ = mat44_inv.topRightCorner<3,1>();
  Eigen::Vector3d pt_b_ = pt_b_r_+pt_b_t_;
//  Eigen::Vector3d pt_b_ = mat44_inv.topLeftCorner<3,3>() * pt_f_ + mat44_inv.topRightCorner<3,1>();

  double diff_mat = (pt-pt_b_).norm();


  std::cout << "diff pt_f: " << (pt_f-pt_f_).norm() << "\n";
  std::cout << "diff pt_b: " << (pt_b-pt_b_).norm() << "\n";
  std::cout << "diff pt_b_r: " << (pt_b_r-pt_b_r_).norm() << "\n";
  std::cout << "diff pt_b_t: " << (pt_b_t-pt_b_t_).norm() << "\n";
  printf("diff: %f, diff_mat: %f\n",diff,diff_mat);
  return 0;
}