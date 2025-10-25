#include <Rcpp.h>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

std::vector<std::pair<int, int>> bresenhamLine(int x1, int y1, int x2, int y2) {
  std::vector<std::pair<int, int>> line;
  int dx = abs(x2 - x1);
  int dy = abs(y2 - y1);
  int sx = (x1 < x2) ? 1 : -1;
  int sy = (y1 < y2) ? 1 : -1;
  int err = dx - dy;
  
  int x = x1;
  int y = y1;
  
  while (true) {
    line.push_back(std::make_pair(x, y));
    
    if (x == x2 && y == y2) break;
    
    int e2 = 2 * err;
    if (e2 > -dy) {
      err -= dy;
      x += sx;
    }
    if (e2 < dx) {
      err += dx;
      y += sy;
    }
  }
  
  return line;
}

// [[Rcpp::export]]
NumericMatrix computeViewshedOptimized(NumericMatrix dsm, NumericMatrix dtm, 
                                       int cx, int cy, double obs_height, double tgt_height) {
  int nrow = dsm.nrow();
  int ncol = dsm.ncol();
  NumericMatrix viewshed(nrow, ncol);
  
  // Observer elevation based on TERRAIN (DTM) + observer height
  double observer_elev = dtm(cx, cy) + obs_height;

#ifdef _OPENMP
  #pragma omp parallel for collapse(2)
#endif
  for (int tx = 0; tx < nrow; ++tx) {
    for (int ty = 0; ty < ncol; ++ty) {
      if (tx == cx && ty == cy) {
        viewshed(tx, ty) = 1;
        continue;
      }

      std::vector<std::pair<int, int>> line = bresenhamLine(cx, cy, tx, ty);
      double max_angle = -1e9;
      bool blocked = false;

      // Check line of sight using SURFACE (DSM) elevations
      for (size_t i = 1; i < line.size() - 1; ++i) {
        int px = line[i].first;
        int py = line[i].second;

        if (px < 0 || px >= nrow || py < 0 || py >= ncol) {
          blocked = true;
          break;
        }

        double dx = px - cx;
        double dy = py - cy;
        double dist = std::sqrt(dx * dx + dy * dy);
        if (dist == 0) continue;

        // Use DSM (surface) elevation for blocking objects
        double angle = (dsm(px, py) - observer_elev) / dist;
        if (angle > max_angle) {
          max_angle = angle;
        }
      }

      if (blocked) {
        viewshed(tx, ty) = 0;
        continue;
      }

      // Target elevation based on TERRAIN (DTM) + target height
      double dx = tx - cx;
      double dy = ty - cy;
      double dist = std::sqrt(dx * dx + dy * dy);
      double target_elev = dtm(tx, ty) + tgt_height;
      double target_angle = (target_elev - observer_elev) / dist;

      viewshed(tx, ty) = (target_angle >= max_angle) ? 1 : 0;
    }
  }

  return viewshed;
}

