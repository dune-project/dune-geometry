// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/geometry/type.hh>

namespace Dune
{

  const unsigned int GeometryType::topologyIdFromVertices[ 5 ][ 16 ]
    = { {  0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u },
        { ~0u,  1u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u },
        { ~0u, ~0u,  0u,  3u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u },
        { ~0u, ~0u, ~0u,  0u,  3u,  4u, ~0u,  7u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u, ~0u },
        { ~0u, ~0u, ~0u, ~0u,  0u,  3u,  4u,  8u,  7u, 11u, ~0u, 12u, ~0u, ~0u, ~0u, 15u } };

} // namespace Dune
