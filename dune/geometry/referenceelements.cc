// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/geometry/referenceelements.hh>

namespace Dune
{

  namespace Impl
  {

    // ReferenceVolumeInverse
    // ----------------------

    unsigned long referenceVolumeInverse ( unsigned int topologyId, int dim )
    {
      assert( (dim >= 0) && (topologyId < numTopologies( dim )) );

      if( dim > 0 )
      {
        unsigned int baseValue = referenceVolumeInverse( baseTopologyId( topologyId, dim ), dim-1 );
        return (isPrism( topologyId, dim ) ? baseValue : baseValue * (unsigned long)dim);
      }
      else
        return 1;
    }

  } // namespace Impl

} // namespace Dune
