#include "../RecoTrueCalculator.hh"

void do_it() {

const std::string respmat_file_name(
    "/uboone/data/users/gardiner/23-sept10-all-universes.root" );

  auto* syst_ptr = new RecoTrueCalculator( respmat_file_name, "../systcalc.conf" );
  auto& syst = *syst_ptr;

  syst.get_covariances();
}
