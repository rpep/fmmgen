#pragma once
#include "utils.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

class vec3 {
  public:
    vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    double norm() {
      return sqrt(norm2());
    }
    double norm2() {
      return x*x + y*y + z*z;
    }

    vec3 operator+(const vec3 &v) {
      return vec3(x + v.x, y + v.y, z + v.z);
    }

    vec3 operator-(const vec3 &v) {
      return vec3(x - v.x, y + v.y, z + v.z);
    }

    vec3& operator+=(const vec3 &v) {
      this->x += v.x;
      this->y += v.y;
      this->z += v.z;
      return *this;
    }

    vec3& operator-=(const vec3 &v) {
      this->x -= v.x;
      this->y -= v.y;
      this->z -= v.z;
      return *this;
    }
  private:
    double x, y, z;
};



/*! \brief Particle class used to store position and dipole moment strength. */
class Particle {
public:
  double x; /*!< x position of the particle */
  double y; /*!< y position of the particle */
  double z; /*!< z position of the particle */
  double q;
  Particle(double x, double y, double z, double q) : x(x), y(y), z(z), q(q) {}
};

class Cell {
public:
  size_t nleaf;  /*!< \brief Number of particles held in cell.

                   This counter is
                   incremented every time a particle is added to it in the
                   \ref build_tree function. This continues to be the case
                   even when the cell has been split, as we use it to keep
                   track of whether a cell has been split or not to
                   save on memory, rather than having another variable. */
  int nchild; /*!< \brief Number of child cells occupied.
                          Binary counter showing whether a given octant is
                          occupied by a child cell.<br>I.e. if 0001001, then
                          there are two child cells held by this cell.
                          DO NOT CHANGE TO size_t!*/
  size_t level; /*!< \brief Level of the tree that the cell sits at.

                  This is 0 for the root cell, 1 for the 1st level, etc.
                  */
  std::vector<size_t> child; /*!< \brief Indices of child octants. */
  std::vector<double> M;
  std::vector<double> L;
  std::vector<size_t> leaf;            /*!< \brief Indices of particles in the cell. */
  double x;            /*!< \brief x coordinates of cell centre. */
  double y;            /*!< \brief y coordinates of cell centre. */
  double z;            /*!< \brief z coordinates of cell centre. */
  double r;            /*!< \brief Radius of cell
                           Must be sufficiently large for the root cell to bound the
                           particles.

                           Note: I may change this in future so it is calculated
                           in build_tree rather than user specified.
                           */
  double rmax;         /*!< \brief maximum distance to particle from cell centre */
  size_t parent; /*!< \brief Index of parent cell of this cell. */
  Cell(double x, double y, double z, double r, size_t parent,
       size_t order, size_t level, size_t ncrit);
  ~Cell();
  Cell(const Cell &other);
  Cell(Cell &&other);
  void clear();
  void resize(size_t order);
  /*! Copy operator for the Cell class */

  Cell &operator=(const Cell &other) {
    this->nleaf = other.nleaf;
    this->nchild = other.nchild;
    this->level = other.level;
    this->child = other.child;
    this->M = other.M;
    this->leaf = other.leaf;
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->r = other.r;
    this->parent = other.parent;
    return *this;
  }
};

void printTreeParticles(std::vector<Cell> &cells, size_t cell,
                        size_t depth);

void add_child(std::vector<Cell> &cells, int octant, size_t p,
               size_t ncrit, size_t order);

void split_cell(std::vector<Cell> &cells, std::vector<Particle> &particles,
                size_t p, size_t ncrit, size_t order);

std::vector<Cell> build_tree(std::vector<Particle> &particles, Cell &root,
                             size_t ncrit, size_t order);
