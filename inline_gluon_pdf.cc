/*! \file
 *  \brief Inline plaquette
 */
#include "chroma.h"
#include "meas/inline/glue/inline_plaquette.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/glue/mesplq.h"
#include "meas/inline/io/named_objmap.h"
#include "inline_gluon_pdf.h"
#include "meas/inline/io/default_gauge_field.h"

#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

#include "inline_gluon_pdf.h"

namespace Chroma 
{  
  //! Plaquette input
  void read(XMLReader& xml, const std::string& path, InlineGluonObsEnv::Params::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);
    param.cgs = CreateGaugeStateEnv::nullXMLGroup();

    switch (version) 
    {
    case 2:
      if (paramtop.count("GaugeState") != 0)
	param.cgs = readXMLGroup(paramtop, "GaugeState", "Name");
      break;

    default:
      QDPIO::cerr << "InlineGluonObsEnv::Params::Param_t: " << version
		  << " unsupported." << std::endl;
      QDP_abort(1);
    }
  }

  //! Plaquette output
  void write(XMLWriter& xml, const std::string& path, const InlineGluonObsEnv::Params::Param_t& param)
  {
    push(xml, path);

    int version = 2;
    write(xml, "version", version);
    xml << param.cgs.xml;

    pop(xml);
  }


  //! Plaquette input
  void read(XMLReader& xml, const std::string& path, InlineGluonObsEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
  }

  //! Plaquette output
  void write(XMLWriter& xml, const std::string& path, const InlineGluonObsEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);

    pop(xml);
  }


  namespace InlineGluonObsEnv
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "GluonOB";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= CreateGaugeStateEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    // Params
    Params::Params() 
    { 
      frequency = 0; 
      param.cgs          = CreateGaugeStateEnv::nullXMLGroup();
      named_obj.gauge_id = InlineDefaultGaugeField::getId();
      xml_file ="";
    }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Params
	read(paramtop, "Param", param);

	// Ids
	read(paramtop, "NamedObject", named_obj);

	// Possible alternate XML file pattern
	if (paramtop.count("xml_file") != 0) {
	  read(paramtop, "xml_file", xml_file);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << "Caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }


    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      if( params.xml_file != "" ) 
      {
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);
	push( xml_out, "Plaquette");
	write(xml_out, "update_no", update_no);
	write(xml_out, "xml_file", xml_file);
	pop(xml_out);

	XMLFileWriter xml(xml_file);
	func(update_no, xml);
      }
      else 
      {
	func(update_no, xml_out);
      }

    }

    void 
    InlineMeas::func(const unsigned long update_no, 
		     XMLWriter& xml_out) 
    {
		
      START_CODE();
    
      // Test and grab a reference to the gauge field
      multi1d<LatticeColorMatrix> u;
      XMLBufferWriter gauge_xml;

      try
      {
	u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);

	// Set the construct state and modify the fields
	{
	  std::istringstream  xml_s(params.param.cgs.xml);
	  XMLReader  gaugetop(xml_s);

	  Handle<CreateGaugeState< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	    cgs(TheCreateGaugeStateFactory::Instance().createObject(params.param.cgs.id,
								    gaugetop,
								    params.param.cgs.path));

	  Handle<GaugeState< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	    state((*cgs)(u));

	  // Pull the u fields back out from the state since they might have been
	  // munged with gaugeBC's
	  u = state->getLinks();
	}
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << InlineGluonObsEnv::name << ": caught dynamic cast error"
		    << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << InlineGluonObsEnv::name << ": std::map call failed: " << e
		    << std::endl;
	QDP_abort(1);
      }

      push(xml_out, "Plaquette");
      write(xml_out, "update_no", update_no);

      Double w_plaq, s_plaq, t_plaq, link; 
      multi2d<Double> plane_plaq;

      MesPlq(u, w_plaq, s_plaq, t_plaq, plane_plaq, link);
      write(xml_out, "w_plaq", w_plaq);
      write(xml_out, "s_plaq", s_plaq);
      write(xml_out, "t_plaq", t_plaq);

      if (Nd >= 2)
      {
	write(xml_out, "plane_01_plaq", plane_plaq[0][1]);
      }

      if (Nd >= 3)
      {
	write(xml_out, "plane_02_plaq", plane_plaq[0][2]);
	write(xml_out, "plane_12_plaq", plane_plaq[1][2]);
      }

      if (Nd >= 4)
      {
	write(xml_out, "plane_03_plaq", plane_plaq[0][3]);
	write(xml_out, "plane_13_plaq", plane_plaq[1][3]);
	write(xml_out, "plane_23_plaq", plane_plaq[2][3]);
      }

      write(xml_out, "link", link);
    
      pop(xml_out); // pop("Plaquette");
    
		multi2d<LatticeColorMatrix> plane_plaq_matrix_1,plane_plaq_matrix_2,plane_plaq_matrix_3,plane_plaq_matrix_4;
		
		plane_plaq_matrix_1.resize(Nd,Nd);
		plane_plaq_matrix_2.resize(Nd,Nd);
		plane_plaq_matrix_3.resize(Nd,Nd);
		plane_plaq_matrix_4.resize(Nd,Nd);

		for(int mu=1; mu < Nd; ++mu)
		{
			for(int nu=0; nu < mu; ++nu)
			{
				plane_plaq_matrix_1[mu][nu]=u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu]);
				plane_plaq_matrix_1[nu][mu]=adj(plane_plaq_matrix_1[mu][nu]);
				plane_plaq_matrix_2[mu][nu]=u[nu]*adj(shift(shift(u[mu],BACKWARD,mu),FORWARD,nu))*adj(shift(u[nu],BACKWARD,mu))*shift(u[mu],BACKWARD,mu);
                                plane_plaq_matrix_2[nu][mu]=adj(plane_plaq_matrix_2[mu][nu]);
                                plane_plaq_matrix_3[mu][nu]=adj(shift(u[mu],BACKWARD,mu))*adj(shift(shift(u[nu],BACKWARD,mu),BACKWARD,nu))*shift(shift(u[mu],BACKWARD,mu),BACKWARD,nu)*shift(u[nu],BACKWARD,nu);
                                plane_plaq_matrix_3[nu][mu]=adj(plane_plaq_matrix_3[mu][nu]);
				plane_plaq_matrix_4[mu][nu]=adj(shift(u[nu],BACKWARD,nu))*shift(u[mu],BACKWARD,nu)*shift(shift(u[nu],BACKWARD,nu),FORWARD,mu)*adj(u[mu]);
                                plane_plaq_matrix_4[nu][mu]=adj(plane_plaq_matrix_4[mu][nu]);
			}
		}
		int lmax=9;
		multi2d<LatticeColorMatrix> FF,FF00,FF01,FF10,FF11;
		FF.resize(Nd,Nd);
		FF00.resize(Nd,Nd);
		FF01.resize(Nd,Nd);
		FF10.resize(Nd,Nd);
		FF11.resize(Nd,Nd);
		for(int mu=1; mu < Nd; ++mu)
		{
			for(int nu=0; nu < mu; ++nu)
			{
				QDPIO::cout << "mu,nu=" << mu << nu << std::endl;
				FF00[mu][nu]=(plane_plaq_matrix_1[mu][nu]-plane_plaq_matrix_1[nu][mu])/2.0;
                                FF01[mu][nu]=(plane_plaq_matrix_2[mu][nu]-plane_plaq_matrix_2[nu][mu])/2.0;
                                FF10[mu][nu]=(plane_plaq_matrix_3[mu][nu]-plane_plaq_matrix_3[nu][mu])/2.0;
                                FF11[mu][nu]=(plane_plaq_matrix_4[mu][nu]-plane_plaq_matrix_4[nu][mu])/2.0;
				
				FF[mu][nu]=(FF00[mu][nu]+FF01[mu][nu]+FF10[mu][nu]+FF11[mu][nu])/4.0;
				FF[nu][mu]=FF[mu][nu];
				FF[mu][mu]=0;
				FF[nu][nu]=0;
			}
		}
		QDPIO::cout << "FF calculated. FF00"<< sum(real(trace(FF00[2][0]*FF00[2][0])))<<", FF01" <<sum(real(trace(FF01[2][0]*FF01[2][0])))<<", FF10"<<sum(real(trace(FF10[2][0]*FF10[2][0])))<<", FF11"<<sum(real(trace(FF11[2][0]*FF11[2][0])))<<", FF"<<sum(real(trace(FF[0][2]*FF[0][2])))<<", FF"<<sum(real(trace(FF[2][0]*FF[2][0]))) << std::endl;	
		multi2d<Double> gluon_o3,gluon_o1,gluon_o2,gluon_o0;
		gluon_o1.resize(Layout::lattSize()[3],lmax);
                gluon_o2.resize(Layout::lattSize()[3],lmax);
                gluon_o3.resize(Layout::lattSize()[3],lmax);
                gluon_o0.resize(Layout::lattSize()[3],lmax);
		multi1d<LatticeColorMatrix> FFz_shift,FFt_shift;
		FFz_shift.resize(Nd);
		FFt_shift.resize(Nd);
		LatticeColorMatrix tmp; 
		for(int mu=0; mu<Nd; ++mu)
			{
			FFz_shift[mu]=FF[mu][2];
			FFt_shift[mu]=FF[mu][3];
			}
		LatticeColorMatrix u_shift;
		QDPIO::cout << "FF assigned"<< std::endl;
		u_shift=u[2];
		LatticeReal op1,op2,op3,op0;
		multi2d<Double> Ftx2, Fty2, Ftz2, Fzx2, Fzy2, FtxFzx, FtyFzy;
		Ftx2.resize(Layout::lattSize()[3],lmax);
		Fty2.resize(Layout::lattSize()[3],lmax);
		Ftz2.resize(Layout::lattSize()[3],lmax);
		Fzx2.resize(Layout::lattSize()[3],lmax);
		Fzy2.resize(Layout::lattSize()[3],lmax);
		FtxFzx.resize(Layout::lattSize()[3],lmax);
		FtyFzy.resize(Layout::lattSize()[3],lmax);
		for(int l=0; l < lmax; ++l)
		{
			QDPIO::cout << "l=" << l << std::endl;
			for(int t=0; t<Layout::lattSize()[3];++t) {gluon_o3[t][l]=0;gluon_o2[t][l]=0;gluon_o1[t][l]=0;gluon_o0[t][l]=0;Ftx2[t][l]=0;Fty2[t][l]=0;Ftz2[t][l]=0;Fzx2[t][l]=0;Fzy2[t][l]=0;FtxFzx[t][l]=0;FtyFzy[t][l]=0;}
			for(int mu=0; mu<Nd; ++mu)
			{
					if(l>0)
					//gluon_o3[l]+=sum(where(Layout::latticeCoordinate(3) == 0,real(trace(FF[mu][2]*u_shift*FF_shift[mu]*adj(u_shift))),LatticeScalar(zero)));
					{op3=real(trace(FF[mu][2]*u_shift*FFz_shift[mu]*adj(u_shift)));
					op0=real(trace(FF[mu][3]*u_shift*FFt_shift[mu]*adj(u_shift)));
					op1=real(trace(FF[mu][2]*u_shift*FFz_shift[mu]*adj(u_shift)));
                                        op2=real(trace(FF[mu][3]*u_shift*FFz_shift[mu]*adj(u_shift)));
					}
					else //gluon_o3[l]+=sum(where(Layout::LatticeCoordinate(3)==0,real(trace(FF[mu][2]*FF[mu][2])),LatticeScalar(zero)));
					{op3=real(trace(FF[mu][2]*FF[mu][2]));
					op0=real(trace(FF[mu][3]*FF[mu][3]));
                                        op1=real(trace(FF[mu][2]*FF[mu][2]));
                                        op2=real(trace(FF[mu][3]*FF[mu][2]));
					}
					multi1d<int> coord;
					coord.resize(4);
		                	for(int i=0; i<Layout::lattSize()[0];++i)
						for(int j=0; j<Layout::lattSize()[1];++j)
							for(int k=0; k<Layout::lattSize()[2];++k)
							for(int t=0; t<Layout::lattSize()[3];++t)
					{coord[0]=i;
					coord[1]=j;
					coord[2]=k;
					coord[3]=t;
					if(mu!=2) gluon_o3[t][l]+=peekSite(op3,coord);
					if(mu!=2&&mu!=3) {
						gluon_o0[t][l]+=peekSite(op0,coord);
						gluon_o1[t][l]+=peekSite(op1,coord);
                                                gluon_o2[t][l]+=peekSite(op2,coord);
					}
if(mu==0){Ftx2[t][l]+=peekSite(op0,coord); Fzx2[t][l]+=peekSite(op1,coord); FtxFzx[t][l]+=peekSite(op2,coord);}
if(mu==1){Fty2[t][l]+=peekSite(op0,coord); Fzy2[t][l]+=peekSite(op1,coord);FtyFzy[t][l]+=peekSite(op2,coord);}
if(mu==2){Ftz2[t][l]+=peekSite(op0,coord);}
					}
					tmp=shift(FFz_shift[mu],FORWARD,2);
					FFz_shift[mu]=tmp;
                                        tmp=shift(FFt_shift[mu],FORWARD,2);
                                        FFt_shift[mu]=tmp;
			QDPIO::cout << gluon_o3[0][l] << std::endl;
			}
			if(l>0) tmp=u[2]*shift(u_shift,FORWARD,2);
			else tmp=u[2];
			u_shift=tmp;
		}
		
		std::string out_fname_c("data_test.dat");
		TextFileWriter fout(out_fname_c);
		fout << "#t z Op1 Op2 Op3 Op4" << "\n";
		for(int l=0; l < lmax; ++l)
		for(int t=0; t<Layout::lattSize()[3];++t)
		{
			fout <<t<<"\t"<<l<<"\t"<< gluon_o0[t][l]<<"\t"<< gluon_o1[t][l]<<"\t"<< gluon_o2[t][l]<<"\t"<< gluon_o3[t][l] << "\n";
		}
		fout.close();
		std::string out_fname_comp("data_test_comp.dat");
                TextFileWriter fout_comp(out_fname_comp);
		fout_comp << "#t z F[tx]^2 F[ty]^2 F[tz]^2 F[xz]^2 F[yz]^2 F[tx]F[zx] F[ty]F[zy]" << "\n";
		for(int l=0; l < lmax; ++l)
                for(int t=0; t<Layout::lattSize()[3];++t)
                {
			fout_comp << t<<"\t"<<l<<"\t"<< Ftx2[t][l]<<"\t"<< Fty2[t][l]<<"\t"<< Ftz2[t][l]<<"\t"<< Fzx2[t][l]<<"\t"<< Fzy2[t][l]<<"\t"<< FtxFzx[t][l]<<"\t"<< FtyFzy[t][l]<<"\n";
		}
		fout_comp.close(); 
      END_CODE();
    } 

  } // namespace InlineGluonObsEnv

} // namespace Chroma
