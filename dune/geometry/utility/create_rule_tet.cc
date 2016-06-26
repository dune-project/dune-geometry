#include <array>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

#include <boost/multiprecision/cpp_dec_float.hpp>

#define CALC_PRECISION  80
#define PRINT_PRECISION 35

using boost::multiprecision::cpp_dec_float;
using mp_type = boost::multiprecision::number<cpp_dec_float<CALC_PRECISION> >;

#define _F(n)  mp_type( #n )
#include "quad-permu.h" // from PHG

#define DIM 3
using BaryCoord = std::array<mp_type, DIM + 1>;
using LocalCoord = std::array<mp_type, DIM>;

// expect k1 <= k2 <= k3 and k1+k2+k3 == p
mp_type exact(int k1, int k2, int k3)
{
  mp_type result = 1;
  int l1 = std::max(1,k1),
      l2 = std::max(1,k2),
      l3 = std::max(1,k3),
      l4 = 3 + k1 + k2 + k3;
  while (std::max(std::max(std::max(l1, l2), l3), l4) > 1)
  {
    result *= mp_type(l1*l2*l3)/mp_type(l4);
    l1 = std::max(1, l1-1);
    l2 = std::max(1, l2-1);
    l3 = std::max(1, l3-1);
    l4 = std::max(1, l4-1);
  }
  return result;
}

// evaluate monomial
template <class T>
T eval(int k1, int k2, int k3, T x, T y, T z)
{
  return pow(x, k1) * pow(y, k2) * pow(z, k3);
}

// calculate residuum of quadrature rules
mp_type residuum(int order,
                 std::vector<mp_type> const& weights,
                 std::vector<LocalCoord> const& points)
{
  mp_type res = 0;
  for (int k1 = 0; k1 < order; ++k1) {
    for (int k2 = k1; k2 < order; ++k2) {
      for (int k3 = k2; k3 < order; ++k3) {
        if (k1 + k2 + k3 > order)
          continue;

        mp_type result = 0;
        for (std::size_t i = 0; i < points.size(); ++i)
          result += weights[i] * eval(k1, k2, k3, points[i][0], points[i][1], points[i][2]);
        res += pow(result - exact(k1, k2, k3), 2);
      }
    }
  }
  return sqrt(res);
}

template <class Out>
void print(Out& out, int m, int order,
           std::vector<mp_type>& weights,
           std::vector<mp_type> const& points_)
{
  std::string reference = "Zhang, Linbo; Cui, Tao; Liu, Hui. A set of symmetric quadrature rules on triangles and tetrahedra. J. Comput. Math. 27 (2009), no. 1, 89--96. See also the referenced library PHG version 0.9.2 that can be found at http://lsec.cc.ac.cn/phg/index_en.htm";

  std::vector<LocalCoord> points;
  points.reserve( points_.size() / (DIM + 1) );

  // convert barycentric coordinates to local coordinates
  for (std::size_t i = 0; i < points_.size(); i += (DIM + 1)) {
    LocalCoord p = {{ points_[i+1], points_[i+2], points_[i+3] }};
    points.push_back(p);
  }
  // weights multiplied with simplex area 1/6
  for (std::size_t i = 0; i < weights.size(); ++i)
    weights[i] /= mp_type(6);

  // add some comments about the quadrature rule
  out << std::setprecision(PRINT_PRECISION) << std::scientific;
  out << "# " << reference << "\n";
  out << DIM << " " << order << "\n";

  // print npoints, points, weights, order in DUNE notation
  for (std::size_t i = 0; i < points.size(); ++i) {
    for (std::size_t j = 0; j < DIM; ++j)
      out << points[i][j] << " ";
    out << weights[i] << "\n";
  }
}


int main(int argc, char** argv)
{
  std::string out_dir = "./rules/tet";
  { // -------------------------------------------------------------------------
    int dim=3, order=2, npoints=4;
    std::vector<mp_type> weights = {
      Dup31(.25)
    };
    std::vector<mp_type> points = {
      Perm31(.13819660112501051517954131656343619)
    };

    std::string fn = out_dir + "/phg_" + std::to_string(order) + "-" + std::to_string(npoints) + ".csv";
    std::ofstream out(fn, std::ios_base::out);
    print(out, npoints, order, weights, points);
  }

  { // -------------------------------------------------------------------------
    int dim=3, order=3, npoints=8;
    std::vector<mp_type> weights = {
      Dup31(.13852796651186214232361769837564129),
      Dup31(.11147203348813785767638230162435871)
    };
    std::vector<mp_type> points = {
      Perm31(.32805469671142664733580581998119743),
      Perm31(.10695227393293068277170204157061650)
    };

    std::string fn = out_dir + "/phg_" + std::to_string(order) + "-" + std::to_string(npoints) + ".csv";
    std::ofstream out(fn, std::ios_base::out);
    print(out, npoints, order, weights, points);
  }

  { // -------------------------------------------------------------------------
    int dim=3, order=4, npoints=14;
    std::vector<mp_type> weights = {
      Dup31(.07349304311636194934358694586367885),
      Dup31(.11268792571801585036501492847638892),
      Dup22(.04254602077708146686093208377328816)
    };
    std::vector<mp_type> points = {
      Perm31(.09273525031089122628655892066032137),
      Perm31(.31088591926330060975814749494040332),
      Perm22(.04550370412564965000000000000000000)
    };

    std::string fn = out_dir + "/phg_" + std::to_string(order) + "-" + std::to_string(npoints) + ".csv";
    std::ofstream out(fn, std::ios_base::out);
    print(out, npoints, order, weights, points);
  }

  { // -------------------------------------------------------------------------
    int dim=3, order=5, npoints=14;
    std::vector<mp_type> weights = {
      Dup31(.11268792571801585079918565233328633),
      Dup31(.07349304311636194954371020548632750),
      Dup22(.04254602077708146643806942812025744)
    };
    std::vector<mp_type> points = {
      Perm31(.31088591926330060979734573376345783),
      Perm31(.09273525031089122640232391373703061),
      Perm22(.04550370412564964949188052627933943)
    };

    std::string fn = out_dir + "/phg_" + std::to_string(order) + "-" + std::to_string(npoints) + ".csv";
    std::ofstream out(fn, std::ios_base::out);
    print(out, npoints, order, weights, points);
  }

  { // -------------------------------------------------------------------------
    int dim=3, order=6, npoints=24;
    std::vector<mp_type> weights = {
      Dup31(.01007721105532064294801323744593686),
      Dup31(.03992275025816749209969062755747998),
      Dup31(.05535718154365472209515327785372602),
      Dup211(.04821428571428571428571428571428571)
    };
    std::vector<mp_type> points = {
      Perm31(.04067395853461135311557944895641006),
      Perm31(.21460287125915202928883921938628499),
      Perm31(.32233789014227551034399447076249213),
      Perm211(.06366100187501752529923552760572698,
              .26967233145831580803409780572760635)
    };

    std::string fn = out_dir + "/phg_" + std::to_string(order) + "-" + std::to_string(npoints) + ".csv";
    std::ofstream out(fn, std::ios_base::out);
    print(out, npoints, order, weights, points);
  }

  { // -------------------------------------------------------------------------
    int dim=3, order=7, npoints=35;
    std::vector<mp_type> weights = {
      Dup4(.09548528946413084886057843611722638),
      Dup31(.04232958120996702907628617079854674),
      Dup22(.03189692783285757993427482408294246),
      Dup211(.03720713072833462136961556119148112),
      Dup211(.00811077082990334156610343349109654)
    };
    std::vector<mp_type> points = {
      Perm4(.25),
      Perm31(.31570114977820279942342999959331149),
      Perm22(.05048982259839636876305382298656247),
      Perm211(.18883383102600104773643110385458576,
              .57517163758700002348324157702230752),
      Perm211(.02126547254148324598883610149981994,
              .81083024109854856111810537984823239)
    };

    std::string fn = out_dir + "/phg_" + std::to_string(order) + "-" + std::to_string(npoints) + ".csv";
    std::ofstream out(fn, std::ios_base::out);
    print(out, npoints, order, weights, points);
  }

  { // -------------------------------------------------------------------------
    int dim=3, order=8, npoints=46;
    std::vector<mp_type> weights = {
      Dup31(.00639714777990232132145142033517302),
      Dup31(.04019044802096617248816115847981783),
      Dup31(.02430797550477032117486910877192260),
      Dup31(.05485889241369744046692412399039144),
      Dup22(.03571961223409918246495096899661762),
      Dup211(.00718319069785253940945110521980376),
      Dup211(.01637218194531911754093813975611913)
    };
    std::vector<mp_type> points = {
      Perm31(.03967542307038990126507132953938949),
      Perm31(.31448780069809631378416056269714830),
      Perm31(.10198669306270330000000000000000000),
      Perm31(.18420369694919151227594641734890918),
      Perm22(.06343628775453989240514123870189827),
      Perm211(.02169016206772800480266248262493018,
              .71993192203946593588943495335273478),
      Perm211(.20448008063679571424133557487274534,
              .58057719012880922417539817139062041)
    };

    std::string fn = out_dir + "/phg_" + std::to_string(order) + "-" + std::to_string(npoints) + ".csv";
    std::ofstream out(fn, std::ios_base::out);
    print(out, npoints, order, weights, points);
  }

  { // -------------------------------------------------------------------------
    int dim=3, order=9, npoints=59;
    std::vector<mp_type> weights = {
      Dup4(.05489853459364812686895885032391298),
      Dup31(.00421825735654367356185795185819147),
      Dup31(.02348412311384798927791501022996111),
      Dup31(.00421283454980389148648831814037819),
      Dup31(.02994712640542812769203037546126163),
      Dup22(.03695441750679136335292416138761121),
      Dup211(.00817349224171051348425319650294732),
      Dup211(.00987978656102278957913113314297149),
      Dup211(.02160718741919244401497646690335203)
    };
    std::vector<mp_type> points = {
      Perm4(.25000000000000000000000000000000000),
      Perm31(.03785502061999503609086515586175707),
      Perm31(.16954439965012220000000000000000000),
      Perm31(.05484140424416689000000000000000000),
      Perm31(.32229717190921058836777748445908171),
      Perm22(.10961777508972033704050355954365052),
      Perm211(.45915766038590539763886410168178216,
              .08004485927247373376034330857923567),
      Perm211(.03296694775357210169727386483414899,
              .71879584022434055051132299796383374),
      Perm211(.18174359672117481549870278661377760,
              .60023700739524674102301240348069459)
    };

    std::string fn = out_dir + "/phg_" + std::to_string(order) + "-" + std::to_string(npoints) + ".csv";
    std::ofstream out(fn, std::ios_base::out);
    print(out, npoints, order, weights, points);
  }

  { // -------------------------------------------------------------------------
    int dim=3, order=10, npoints=79;
    std::vector<mp_type> weights = {
      Dup4(.04574189830483037077884770618329337),
      Dup31(.01092727610912416907498417206565671),
      Dup31(.00055352334192264689534558564012282),
      Dup31(.02569337913913269580782688316792080),
      Dup22(.00055387649657283109312967562590035),
      Dup211(.01044842402938294329072628200105773),
      Dup211(.02513844602651287118280517785487423),
      Dup211(.01178620679249594711782155323755017),
      Dup211(.01332022473886650471019828463616468),
      Dup211(.00615987577565961666092767531756180)
    };
    std::vector<mp_type> points = {
      Perm4(.25000000000000000000000000000000000),
      Perm31(.11425191803006935688146412277598412),
      Perm31(.01063790234539248531264164411274776),
      Perm31(.31274070833535645859816704980806110),
      Perm22(.01631296303281644000000000000000000),
      Perm211(.03430622963180452385835196582344460,
              .59830121060139461905983787517050400),
      Perm211(.12346418534551115945916818783743644,
              .47120066204746310257913700590727081),
      Perm211(.40991962933181117418479812480531207,
              .16546413290740130923509687990363569),
      Perm211(.17397243903011716743177479785668929,
              .62916375300275643773181882027844514),
      Perm211(.03002157005631784150255786784038011,
              .81213056814351208262160080755918730)
    };

    std::string fn = out_dir + "/phg_" + std::to_string(order) + "-" + std::to_string(npoints) + ".csv";
    std::ofstream out(fn, std::ios_base::out);
    print(out, npoints, order, weights, points);
  }

  { // -------------------------------------------------------------------------
    int dim=3, order=11, npoints=96;
    std::vector<mp_type> weights = {
      Dup31(.01612698613577620369120244222737879),
      Dup31(.00178872341812357138976990346996962),
      Dup31(.00847529348343123401863799968389086),
      Dup31(.01238021263944669050859562763135516),
      Dup31(.02205586697199415746140963638568037),
      Dup31(.02295765467664274421265594265203307),
      Dup22(.00120553827014535727045055662252294),
      Dup22(.02479381575164443454447803302296997),
      Dup211(.01203878836480353606935457416590660),
      Dup211(.00189370204498242146248858917618493),
      Dup211(.01838752922255814184581020943433469),
      Dup211(.00375249249801662461193260176157591),
      Dup211(.00633289841693951300885921328914879)
    };
    std::vector<mp_type> points = {
      Perm31(.12460560449278830000000000000000000),
      Perm31(.02609630765687464746851542316261877),
      Perm31(.07193883255798884087330011042809557),
      Perm31(.32611122454203676937273102302894204),
      Perm31(.29405882789858127213310307732130217),
      Perm31(.19271399104965490000000000000000000),
      Perm22(.00047127204692773946587837159205225),
      Perm22(.10321360207480949336085123341390539),
      Perm211(.04349989920159741251267172033621503,
              .63045319723555591476353398203997141),
      Perm211(.01414839289422299290755441603794058,
              .82491678632147090000000000000000000),
      Perm211(.21646077368258425486341884576246642,
              .52711130286496480000000000000000000),
      Perm211(.13301884366834711587538262083530116,
              .73318551371398651551736762818473584),
      Perm211(.44054756810613723082959230959880706,
              .11506799584377921703650823955291194)
    };

    std::string fn = out_dir + "/phg_" + std::to_string(order) + "-" + std::to_string(npoints) + ".csv";
    std::ofstream out(fn, std::ios_base::out);
    print(out, npoints, order, weights, points);
  }

  { // -------------------------------------------------------------------------
    int dim=3, order=12, npoints=127;
    std::vector<mp_type> weights = {
      Dup4(.02340581914868067999082580773836836),
      Dup31(.00484469946470415656870798306091558),
      Dup31(.00079865303812732982185563521014343),
      Dup31(.01311872008808756207964488505025527),
      Dup31(.02352182961292765917274505054313770),
      Dup31(.00210860882494149803857437048649497),
      Dup31(.00047839298963616600187228601742259),
      Dup22(.00204546234216855322941711800170502),
      Dup211(.00334576331671817115245418532677178),
      Dup211(.01181044822479275264785338274950585),
      Dup211(.00290156990282342152841364375092118),
      Dup211(.00949250645501753676094846901252898),
      Dup211(.02094018358085748583183796760479700),
      Dup211(.00171435866337409051521874943702732),
      Dup1111(.00759915954173370886076474450830409)
    };
    std::vector<mp_type> points = {
      Perm4(.25000000000000000000000000000000000),
      Perm31(.19318721110347230000000000000000000),
      Perm31(.01811701371436566878506928822499717),
      Perm31(.10700751831426066518406159227423033),
      Perm31(.29936173715970702940603127680004538),
      Perm31(.33333033333333333042835213613025030),
      Perm31(.16575369007421640000000000000000000),
      Perm22(.04009986052352575650366980228640728),
      Perm211(.01951844463761131301132122485607343,
              .59982639757597731668263005976738196),
      Perm211(.24970741896308715787490891769354198,
              .47400425629911050000000000000000000),
      Perm211(.07674205857869954726322831328843659,
              .83056291375422969598432041821082569),
      Perm211(.43011409627915217536723647418133112,
              .02265922072588833582931396831630072),
      Perm211(.12197854304894211937147375564906792,
              .47765370899783134571567376444973682),
      Perm211(.01480482319031682427540691439704854,
              .81083799468092699988474915243749073),
      Perm1111(.65250697573013212016385330106711095,
               .22646235632397177636617160407210034,
               .02251830769546778956654013747639605)
    };

    std::string fn = out_dir + "/phg_" + std::to_string(order) + "-" + std::to_string(npoints) + ".csv";
    std::ofstream out(fn, std::ios_base::out);
    print(out, npoints, order, weights, points);
  }

}
