// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_REFINEMENT_VIRTUALREFINEMENTTAGS_HH
#define DUNE_GEOMETRY_REFINEMENT_VIRTUALREFINEMENTTAGS_HH

namespace Dune{
  namespace VirtualRefinementTag{

    class Intervals{
      int intervals_=1;
      // forbid implicit conversion from int
    public:
      explicit Intervals(int i) : intervals_(i) {}

      int intervals() const { return intervals_; }
    };

    class Levels: public Intervals{
    public:
      explicit Levels(int l) : Intervals(1<<l) {}
    };

  }
}

#endif //DUNE_GEOMETRY_REFINEMENT_VIRTUALREFINEMENTTAGS_HH
