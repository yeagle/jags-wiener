#include <Module.h>
#include <distributions/DWiener.h>

using std::vector;

namespace wiener {

    class WIENERModule : public Module {
    public:
	WIENERModule();
	~WIENERModule();
    };

    WIENERModule::WIENERModule() 
	: Module("wiener")
    {

	//Load distributions
	insert(new DWiener);

  //Load functions 
  /*
	insert(new DFunction(dist));
	insert(new PFunction(dist));
	insert(new QFunction(dist));
  */

    }

    WIENERModule::~WIENERModule() {

	vector<Distribution*> const &dvec = distributions();
	for (unsigned int i = 0; i < dvec.size(); ++i) {
	    delete dvec[i];
	}
    }

}

wiener::WIENERModule _wiener_module;
