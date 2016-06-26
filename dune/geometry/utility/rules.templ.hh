// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_QUADRATURE_${name | uppercase}_HH
#define DUNE_GEOMETRY_QUADRATURE_${name | uppercase}_HH

#ifndef DUNE_INCLUDING_IMPLEMENTATION
#error This is a private header that should not be included directly.
#error Use #include <dune/geometry/quadraturerules.hh> instead.
#endif

namespace Dune {

  /************************************************
   * Quadraturerule for ${geometry}
   *************************************************/

  /** \brief Quadrature rules for ${geometry}
      \ingroup Quadrature
   */
  template<typename ct, int dim>
  class ${name}QuadratureRule;

  /** \brief Quadrature rules for ${geometry}
      \ingroup Quadrature
   */
  template<typename ct>
  class ${name}QuadratureRule<ct, ${dim}> : public QuadratureRule<ct, ${dim}>
  {
  public:
    /** \brief The highest quadrature order available */
    enum { highest_order = ${maxorder} };
  private:
    friend class QuadratureRuleFactory<ct, ${dim}>;
    ${name}QuadratureRule (int p);
    ~${name}QuadratureRule() {}
  };

  template<typename ct>
  ${name}QuadratureRule<ct, ${dim}>::${name}QuadratureRule(int order)
    : QuadratureRule<ct, ${dim}>(GeometryType(GeometryType::${geometry}, ${dim}))
  {
    switch (order)
    {
    %for o in range(0, maxorder+1):
      case ${o}:
<%    rule = rules[o] %>\
      %if rule and rule['points'] and len(rule['points']) > 0:
<%      points, weights = rule['points'], rule['weights'] %>\
        // Symmetric rule: dim = ${dim}, order = ${o}, npoints = ${len(points)}
        % if len(rule['reference']) > 0:
        // Source:
        // ${("\n" + " "*8 + "// ").join(rule['reference'])}
        % endif
        this->delivered_order = ${o};
        this->reserve(${len(points)});
        this->insert(this->end(), {
        %for i,p in enumerate(points):
          {
            { ${(",\n" + " "*14).join(map(cast, map(str, p)))} },
            ${str(weights[i]) | cast}
          ${"}," if i < len(points)-1 else "}"}
        %endfor
        });
        break;

      %endif
    %endfor
      default:
        DUNE_THROW(QuadratureOrderOutOfRange,
                  "QuadratureRule for order " << order << " and GeometryType "
                                              << this->type() << " not available");
    }
  }

} // end namespace Dune

#endif // DUNE_GEOMETRY_QUADRATURE_${name | uppercase}_HH
