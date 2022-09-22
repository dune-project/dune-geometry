// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <dune/geometry/referenceelementimplementation.hh>

namespace Dune
{

  namespace Geo
  {

    namespace Impl
    {

      // size
      // ----

      unsigned int size ( unsigned int topologyId, int dim, int codim )
      {
        assert( (dim >= 0) && (topologyId < numTopologies( dim )) );
        assert( (0 <= codim) && (codim <= dim) );

        if( codim > 0 )
          {
            const unsigned int baseId = baseTopologyId( topologyId, dim );
            const unsigned int m = size( baseId, dim-1, codim-1 );

            if( isPrism( topologyId, dim ) )
              {
                const unsigned int n = (codim < dim ? size( baseId, dim-1, codim ) : 0);
                return n + 2*m;
              }
            else
              {
                assert( isPyramid( topologyId, dim ) );
                const unsigned int n = (codim < dim ? size( baseId, dim-1, codim ) : 1);
                return m+n;
              }
          }
        else
          return 1;
      }



      // subTopologyId
      // -------------

      unsigned int subTopologyId ( unsigned int topologyId, int dim, int codim, unsigned int i )
      {
        assert( i < size( topologyId, dim, codim ) );
        const int mydim = dim - codim;

        if( codim > 0 )
          {
            const unsigned int baseId = baseTopologyId( topologyId, dim );
            const unsigned int m = size( baseId, dim-1, codim-1 );

            if( isPrism( topologyId, dim ) )
              {
                const unsigned int n = (codim < dim ? size( baseId, dim-1, codim ) : 0);
                if( i < n )
                  return subTopologyId( baseId, dim-1, codim, i ) | ((unsigned int)prismConstruction << (mydim - 1));
                else
                  return subTopologyId( baseId, dim-1, codim-1, (i < n+m ? i-n : i-(n+m)) );
              }
            else
              {
                assert( isPyramid( topologyId, dim ) );
                if( i < m )
                  return subTopologyId( baseId, dim-1, codim-1, i );
                else if( codim < dim )
                  return subTopologyId( baseId, dim-1, codim, i-m ) | ((unsigned int)pyramidConstruction << (mydim - 1));
                else
                  return 0u;
              }
          }
        else
          return topologyId;
      }



      // subTopologyNumbering
      // --------------------

      void subTopologyNumbering ( unsigned int topologyId, int dim, int codim, unsigned int i, int subcodim,
                                  unsigned int *beginOut, unsigned int *endOut )
      {
        assert( (codim >= 0) && (subcodim >= 0) && (codim + subcodim <= dim) );
        assert( i < size( topologyId, dim, codim ) );
        assert( (endOut - beginOut) == size( subTopologyId( topologyId, dim, codim, i ), dim-codim, subcodim ) );

        if( codim == 0 )
          {
            for( unsigned int j = 0; (beginOut + j) != endOut; ++j )
              *(beginOut + j) = j;
          }
        else if( subcodim == 0 )
          {
            assert( endOut = beginOut + 1 );
            *beginOut = i;
          }
        else
          {
            const unsigned int baseId = baseTopologyId( topologyId, dim );

            const unsigned int m = size( baseId, dim-1, codim-1 );

            const unsigned int mb = size( baseId, dim-1, codim+subcodim-1 );
            const unsigned int nb = (codim + subcodim < dim ? size( baseId, dim-1, codim+subcodim ) : 0);

            if( isPrism( topologyId, dim ) )
              {
                const unsigned int n = size( baseId, dim-1, codim );
                if( i < n )
                  {
                    const unsigned int subId = subTopologyId( baseId, dim-1, codim, i );

                    unsigned int *beginBase = beginOut;
                    if( codim + subcodim < dim )
                      {
                        beginBase = beginOut + size( subId, dim-codim-1, subcodim );
                        subTopologyNumbering( baseId, dim-1, codim, i, subcodim, beginOut, beginBase );
                      }

                    const unsigned int ms = size( subId, dim-codim-1, subcodim-1 );
                    subTopologyNumbering( baseId, dim-1, codim, i, subcodim-1, beginBase, beginBase+ms );
                    for( unsigned int j = 0; j < ms; ++j )
                      {
                        *(beginBase+j) += nb;
                        *(beginBase+j+ms) = *(beginBase+j) + mb;
                      }
                  }
                else
                  {
                    const unsigned int s = (i < n+m ? 0 : 1);
                    subTopologyNumbering( baseId, dim-1, codim-1, i-(n+s*m), subcodim, beginOut, endOut );
                    for( unsigned int *it = beginOut; it != endOut; ++it )
                      *it += nb + s*mb;
                  }
              }
            else
              {
                assert( isPyramid( topologyId, dim ) );

                if( i < m )
                  subTopologyNumbering( baseId, dim-1, codim-1, i, subcodim, beginOut, endOut );
                else
                  {
                    const unsigned int subId = subTopologyId( baseId, dim-1, codim, i-m );
                    const unsigned int ms = size( subId, dim-codim-1, subcodim-1 );

                    subTopologyNumbering( baseId, dim-1, codim, i-m, subcodim-1, beginOut, beginOut+ms );
                    if( codim+subcodim < dim )
                      {
                        subTopologyNumbering( baseId, dim-1, codim, i-m, subcodim, beginOut+ms, endOut );
                        for( unsigned int *it = beginOut + ms; it != endOut; ++it )
                          *it += mb;
                      }
                    else
                      *(beginOut + ms) = mb;
                  }
              }
          }
      }



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

  } // namespace Geo

} // namespace Dune
