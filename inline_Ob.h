#ifndef __inline_Ob_h__
#define __inline_Ob_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"


namespace Chroma 
{
  // A namespace for this particular   measurement
  namespace InlineObEnv 
  {
    extern const std::string name;
    bool registerAll();
  

    //! Parameter structure
    struct InlineObParams 
    {
	// Default constructor -- forward declaration
	InlineObParams();
	
	// Construct from XML -- forward declaration
	InlineObParams(XMLReader& xml_in, const std::string& path);
	
	// Write myself out
	void write(XMLWriter& xml_out, const std::string& path);

	struct NamedObject_t
	{
	    std::string gauge_id;
	} named_obj;

	struct Src_t
	{
	    multi1d<int> srcLoc;
	    int t_start;
	    int t_end;
	};
	
	multi1d<Src_t> srcs;
	int radius;
	double frequency;

	int max_mom2 ; // max p^2

	
	// Parameters to pull from the xml file
	std::string xml_file; //optional xml output

		
    }; // struct
  }; // namespace InlineObEnv

  class InlineMyMeas : public AbsInlineMeasurement 
  {
  public:
      // Constructor -- default -- do nothing
      ~InlineMyMeas() {}
      
      // Constructor -- from param struct -- copy param struct inside me
      InlineMyMeas(const InlineObEnv::InlineObParams& p) : params(p) {}
      
      // Constructor -- copy constructor -- copy param struct from argument
      InlineMyMeas(const InlineMyMeas& p) : params(p.params) {}

    // Boiler plate
    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 


  private:
    void func(const unsigned long update_no,
	      XMLWriter& xml_out);

    //function to get the time direction
    static int tDir() { return Nd-1;}

    //Get if the location is valid given Radius
    bool validLocation(const multi1d<int>& t_coords,
		       const multi1d<int>& t_src,
		       int R);

    //Returns unnormalized O_b at coordinates coords
    //Also returns E^2 and B^2 at that point
    Double getO_b(const multi1d<int>& coords,
		  const multi2d<LatticeColorMatrix>& plane_plaq,
		  Double &E, Double &B);
    
    //Vector contains O_b found at each time t for passed src
    // Assumes tDir == 3
    void getO_b(std::vector<Double>& O_b,
		const  InlineObEnv::InlineObParams::Src_t src,
		const int radius,
		const LatticeColorMatrix& Oa0);
    
    ColorMatrix get_G(const multi1d<int>& coords, int mu, int nu,
		      const multi2d<LatticeColorMatrix>& P);

    int leviCh(int i, int j, int k)
    {
      if(i == j || i == k || j == k)
	return 0;
      if( (i+1)%3 == j%3)
	return 1;
      return -1;
    };
    

    //private parameters
    InlineObEnv::InlineObParams params;
    
  };

};

namespace QDP
{
    inline void read(XMLReader& xml, const std::string& path, Chroma::InlineObEnv::InlineObParams::Src_t& src)
    {
	XMLReader namedTop(xml, path);

	read(namedTop, "t_src", src.srcLoc);
	read(namedTop, "t_start", src.t_start);
	read(namedTop, "t_end", src.t_end);
	
    }

    inline void write(XMLWriter& xml, const std::string& path, const Chroma::InlineObEnv::InlineObParams::Src_t& src)
    {
	push(xml, path);

	write(xml, "t_src", src.srcLoc);
	write(xml, "t_start", src.t_start);
	write(xml, "t_end", src.t_end);
	
	pop(xml);
    }
    inline void write(QDP::XMLWriter& xml, const std::string& path, const Chroma::InlineObEnv::InlineObParams::NamedObject_t& named_obj)
    {
	push(xml, path);
	write(xml, "gauge_id", named_obj.gauge_id);
	pop(xml);
    }

    /** Functions to read input xml file **/
    inline void read(XMLReader& xml, const std::string& path, Chroma::InlineObEnv::InlineObParams::NamedObject_t& named_obj)
    {
	XMLReader namedTop(xml, path);
	
	read(namedTop, "gauge_id", named_obj.gauge_id);
    }

}

#endif
