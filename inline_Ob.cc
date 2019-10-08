#include "inline_Ob.h"
#include "chroma.h"
#include <math.h>
#include <complex>
#include <iostream>
#include <string>
#include "stdio.h"



namespace Chroma 
{
    
    namespace InlineObEnv 
    {
	//Name of the measurement to be called in the XML input file
	const std::string name = "GMF_O_b";
	
	//This function is used with the factory thing
	AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
						const std::string& path) 
	{
	    //Create new instance of measurement class using params
	    //from the passed xml file
	    return new InlineMyMeas(InlineObParams(xml_in, path));
	}
		
	// Local registration flag
	namespace {
	    bool registered = false;
	}
	
	// Function to register all the factories
	bool registerAll() 
	{
	    bool success = true; 
	    if (! registered)
	    {
		success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
		QDPIO::cout << "Registering " << name << " " << success << std::endl;
		registered = true;
	    }
	    return success;
	}
	
	
	
	/*** Implementation of Parameter functions ***/
	
	//Set default parameters
	InlineObParams::InlineObParams() { frequency = 0; radius = 0; srcs.resize(1); }
	
	//Read parameters in from xml file
	InlineObParams::InlineObParams(XMLReader& xml_in, const std::string& path) 
	{
	    try 
	    {
		XMLReader paramtop(xml_in, path);
		
		if (paramtop.count("Frequency") == 1)
		    read(paramtop, "Frequency", frequency);
		else
		    frequency = 1;
		
		read(paramtop, "NamedObject", named_obj);

		//Read in the starting position time (ie the time
		//of the first source) and the time between each
		//source. This assumes that each source is equally
		//spaced in time. See branch First-O_b for code
		//that drops this assumtion and takes sources
		//as the arguments (loc & start/stop times)
		read(paramtop, "Multi_Src", srcs);

		if(paramtop.count("radius") == 1)
		    read(paramtop, "radius", radius);
		else
		    radius = 0;

		
		// Possible alternate XML file pattern
		if (paramtop.count("xml_file") != 0) 
		{
		    read(paramtop, "xml_file", xml_file);
		} else
		    xml_file = "";
		
	    }
	    catch(const std::string& e) 
	    {
		QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
		QDP_abort(1);
	    }
	}
	
	
	// Write loaded params to xml file
	void InlineObParams::write(XMLWriter& xml_out, const std::string& path) 
	{
	    push(xml_out, path);
	    
	    // write our all params
	    QDP::write(xml_out, "Multi_Src", srcs);
	    //QDP::write(xml_out, "Named_Object", named_obj);
	    QDP::write(xml_out, "radius", radius);

		  
	    
	    if(xml_file != "")
		QDP::write(xml_out, "xml_file", xml_file);
	    pop(xml_out);
	}
    } //End namespace InlineObEnv

    /*** Inline Measurement function implimentation ***/
    // Function call
    void InlineMyMeas::operator()(unsigned long update_no,
				  XMLWriter& xml_out) 
    {
	// If xml file not empty, then use alternate
	if (params.xml_file != "")
	{
	    std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	    push(xml_out, "GMF_O_b");
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
    
    
    /*** Measurement code stars here ***/
    void InlineMyMeas::func(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
	START_CODE();
	
	QDPIO::cout << InlineObEnv::name << ": Begining" << std::endl;
	
	StopWatch snoop;
	snoop.reset();
	snoop.start();

	//Print out boilerplate stuff to xml
	push(xml_out, "GMF_O_b");
	write(xml_out, "update_no", update_no);

	//Write out the input
	params.write(xml_out, "Input");


	/** Calculate the two dimensional plaquettes **/
	
	//Get link matrix
	multi1d<LatticeColorMatrix> u;
	u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

	//Create variables to store the plaquette and the average
	//trace of the plaquette
	multi3d<LatticeColorMatrix> plane_plaq;
	multi2d<Double> tr_plane_plaq;
	plane_plaq.resize(4,Nd,Nd); //Quadrent, plane
	tr_plane_plaq.resize(Nd,Nd);
	Double w_plaq;
	Double s_plaq;
	Double t_plaq;

	/* Calculate the plaquette and the average trace of
	 * the plaquette
	 */
	for(int mu = 0; mu < Nd; mu++)
	{
	    for(int nu = mu+1; nu < Nd; nu++)
	    {
		LatticeColorMatrix tmp, tmp2, tmp3;
		//Do first quadrent
		tmp = shift(u[nu], FORWARD, mu);
		tmp2 = u[mu] * tmp;
		tmp = shift(u[mu], FORWARD, nu);
		tmp3 = u[nu] * tmp;

		//Record the plaquette in a cross section
		plane_plaq[0][nu][mu] = tmp2*adj(tmp3);

		//Do second quadrent
		tmp = shift(shift(u[mu],FORWARD, nu), BACKWARD, mu);
		tmp2 = u[nu] * adj(tmp);
		tmp = shift(u[nu], BACKWARD, mu);
		tmp3 = adj(tmp) * shift(u[mu], BACKWARD, mu);

		plane_plaq[1][nu][mu] = tmp2*tmp3;

		//Do third quad
		tmp = shift(u[mu], BACKWARD, mu);
		//tmp2 = shift(shift(u[nu], BACKWARD, nu), BACKWARD, mu);
		//tmp2 = adj(tmp) * adj(tmp2);
                tmp2 = adj(tmp)*adj(shift(shift(u[nu], BACKWARD, nu), BACKWARD, mu));
		tmp = shift(shift(u[mu], BACKWARD, nu), BACKWARD, mu);
		tmp3 = tmp * shift(u[nu], BACKWARD, nu);

		plane_plaq[2][nu][mu] = tmp2*tmp3;

		//Do fourth quad
		tmp = shift(u[nu], BACKWARD, nu);
		tmp2 = adj(tmp) * shift(u[mu], BACKWARD, nu);
		tmp3 = shift(shift(u[nu], FORWARD, mu), BACKWARD, nu) * adj(u[mu]);

		plane_plaq[3][nu][mu] = tmp2*tmp3;
		
		tr_plane_plaq[nu][mu] = sum(real(trace(plane_plaq[0][nu][mu])));
		
		//Normalize the plane
		tr_plane_plaq[nu][mu] /= Double(Layout::vol() * Nc);
		plane_plaq[0][mu][nu] = plane_plaq[0][nu][mu]; //symmetric
		tr_plane_plaq[mu][nu] = tr_plane_plaq[nu][mu]; //symmetric
		
		w_plaq += tr_plane_plaq[nu][mu];

		//record the time/space plaq
		if(nu == tDir())
		    t_plaq += tr_plane_plaq[nu][mu];
		else
		    s_plaq += tr_plane_plaq[nu][mu];

	    } //end nu loop
	} //Found normalized plane and unorm plaq: end mu loop

	//Normalize the average plaquette for space/time/whole
	w_plaq *= 2 / Double(Nd*(Nd-1));
	t_plaq /= (Nd - 1);
	if(Nd > 2)
	    s_plaq *= 2 / Double((Nd-2)*(Nd-1));
	
        QDPIO::cout << w_plaq << " " << s_plaq << " " << t_plaq << std::endl;

	/** Record the results in the xml file **/
	write(xml_out, "w_plaq", w_plaq);
	write(xml_out, "s_plaq", s_plaq);
	write(xml_out, "t_plaq", t_plaq);

	/** Write plane plaq to xml file **/
	
	for(int mu = 0; mu < Nd; mu++)
	{
	    for(int nu = mu+1; nu < Nd; nu++)
	    {
		write(xml_out, "plane_plaq_" + std::to_string(mu) +
		      std::to_string(nu), tr_plane_plaq[mu][nu]);
	    }
	}
	

        LatticeColorMatrix OA;
        for(int z = 0; z < 3; z++)
        {
	QDPIO::cout << "Finding F" << std::endl;
        int nn = 8;
	/** Find F_{n,munu} **/
	multi2d<LatticeColorMatrix> F, Fn, nF;
//        multi3d<LatticeColorMatrix> Fi;
	F.resize(Nd,Nd);
        Fn.resize(Nd,Nd);
        nF.resize(Nd,Nd);

//        nF.resize(Nd,Nd);
//        Fi.resize(4,Nd,Nd);
//        int z = 3;
        LatticeColorMatrix u0, un, u1, nu;
        u0 = u[z];
        un = u[z];
        u1 = adj(shift(u[z], BACKWARD, z));
        nu = u1;

        for(int mu = 0; mu < Nd; mu++)
        {
            for(int nu = mu+1; nu < Nd; nu++)
            {
                for(int i = 0; i < 4; i++)
                {
                    F[nu][mu] += plane_plaq[i][nu][mu] - adj(plane_plaq[i][nu][mu]);
                }
                F[nu][mu] *= 1/8.0;
                F[mu][nu]= -F[nu][mu];
                Fn[nu][mu] = F[nu][mu];
                Fn[mu][nu] = -Fn[nu][mu];
                nF[nu][mu] = F[nu][mu];
                nF[mu][nu] = -nF[nu][mu];
            }
        }    


        for(int n = 1; n < nn; n++)
        {
	for(int mu = 0; mu < Nd; mu++)
	{
	    for(int nu = mu+1; nu < Nd; nu++)
	    {
//                nF[nu][mu] = F[nu][mu];
//                Fni[0][nu][mu]=F[nu][mu];
//                Fni[0][mu][nu]=Fni[0][nu][mu];
                Fn[nu][mu] = shift(Fn[nu][mu], FORWARD, z);
                Fn[mu][nu] = -Fn[nu][mu];
                nF[nu][mu] = shift(nF[nu][mu], BACKWARD, z);
                nF[mu][nu] = -nF[nu][mu];

//                       Fni[j][nu][mu]= Fn[nu][mu];
//                       Fni[j][mu][nu]= Fni[j][nu][mu];
//                       nF[nu][mu] = shift(nF[nu][mu], BACKWARD, z);
//                       Fni[nu][mu][j+nn-1]= nF[nu][mu];
//                       Fni[mu][nu][j+nn-1]= Fni[nu][mu][j+nn-1];
            }
        }
              
        if(n!=1)
        {
          u0 =  shift(u0, FORWARD, z);
          un = un*u0;
          u1 =  shift(u1, BACKWARD, z);
          nu = nu*u1;

        }

        QDPIO::cout << "Finding E/B" << std::endl;


        multi1d<LatticeColorMatrix> E,B,En,Bn,nE,nB;

        E.resize(3);
        B.resize(3);
        En.resize(3);
        Bn.resize(3);
        nE.resize(3);
        nB.resize(3);


        LatticeColorMatrix O0, O00, O1, O10, O2, O20, O3, O30, Oa, Oa0, Od, Od0, OC, OC0, OD, OD0;

        O0 *= 0.0;
        O00 *= 0.0;
        O1 *= 0.0;
        O10 *= 0.0;
        O2 *= 0.0;
        O20 *= 0.0;
        O3 *= 0.0;
        O30 *= 0.0;
        Oa *= 0.0;
        Oa0 *= 0.0;
        Od *= 0.0;
        Od0 *= 0.0;
        OC *= 0.0;
        OC0 *= 0.0;
        OD *= 0.0;
        OD0 *= 0.0;
        
        //manually do the z=2 operators
        u2 = shift(u[z], FORWARD, z);
        un = u[z]*u2;
        for(int mu = 0; mu < Nd; mu++)
        {
            for(int nu = mu+1; nu < Nd; nu++)
            {
                Fn[nu][mu] = shift(shift(F[nu][mu], FORWARD, z), FORWARD, z);
                Fn[mu][nu] = -Fn[nu][mu];
            }
        }
        //manually do the z=1 operators, and check the shift function
/*
        u2 = shift(shift(u[z], FORWARD, z), BACKWARD, z);
        un = u2;
        for(int mu = 0; mu < Nd; mu++)
        {
            for(int nu = mu+1; nu < Nd; nu++)
            {
                Fn[nu][mu] = shift(shift(shift(F[nu][mu], FORWARD, z), FORWARD, z), BACKWARD,z);
                Fn[mu][nu] = -Fn[nu][mu];
            }
        }
*/

		
        for(int i=0;i<4; i++)
        {
            O3 += F[z][i]*un*Fn[z][i]*adj(un);
            //O3 += F[z][i]*un*Fn[z][i]*adj(u2)*adj(u[z]);
            O30 += F[z][i]*F[z][i];
        }




        for(int i=0;i<4; i++)
        {
            O0 += F[3][i]*un*Fn[3][i]*adj(un);
            O00 += F[3][i]*F[3][i];
        }


        for(int i=0;i<4; i++)
        {
           for(int j=0;j<i;j++)
              {
                   O0 -= 0.5*F[j][i]*un*Fn[j][i]*adj(un);
                   O00 -= 0.5*F[j][i]*F[j][i];
              }
        }



        for(int i=0;i<4; i++)
        {
            O2 += F[z][i]*un*Fn[z][i]*adj(un);
            O20 += F[z][i]*F[z][i];
        }


        for(int i=0;i<4; i++)
        {
           for(int j=0;j<i;j++)
              {
                   O2 -= 0.5*F[j][i]*un*Fn[j][i]*adj(un);
                   O20 -= 0.5*F[j][i]*F[j][i];
              }
        }


        for(int i=0;i<4; i++)
        {
              O1 += F[3][i]*un*Fn[z][i]*adj(un);
              O10 += F[3][i]*F[z][i];
        }


        for(int i=0;i<4; i++)
        {
              Oa += F[3][i]*un*Fn[3][i]*adj(un);
              Oa0 += F[3][i]*F[3][i];
        }

        for(int i=0;i<3; i++)
        {
            if(i!=z)
            {
               Od += F[z][i]*un*Fn[z][i]*adj(un);
               Od0 += F[z][i]*F[z][i];
            } 
        }

        for(int i = 0; i < 3; i++)
          {
            E[i] = F[3][i];
            En[i] = Fn[3][i];
            nE[i] = nF[3][i];
            int j=(i+1)%3;
            int k=(i+2)%3;
            B[i] += leviCh(i,j,k)*F[j][k];
            Bn[i] += leviCh(i,j,k)*Fn[j][k];
            nB[i] += leviCh(i,j,k)*nF[j][k];
          }
        std::vector<Double> VecO0, VecO00, VecO1, VecO10, VecO2, VecO20, VecO3, VecO30, VecOa, VecOa0, VecOd, VecOd0, VecOC, VecOC0, VecOD, VecOD0;
        QDPIO::cout << "Finding Oa" << std::endl;
//        LatticeColorMatrix Oa, Oam, Oa0, Oc, Oc0, Od, Od0;
        for(int i = 0; i < 3;i++)
        {   //O00 += E[i]*E[i]- B[i]*B[i];
            //O0 += E[i]*un*En[i]*adj(un)-B[i]*un*Bn[i]*adj(un);
//            Oam += E[i]*nu*nE[i]*adj(nu)-B[i]*nu*nB[i]*adj(nu);

        
            if(i!=z)
            {
            OC0 += E[i]*B[i];
            OC += E[i]*un*Bn[i]*adj(un);
            OD0 += E[i]*nu*nB[i]*adj(nu);
            OD += B[i]*un*En[i]*adj(un);

            }

        }

        QDPIO::cout << "Finding Tr(Oa)" << std::endl;
        for(int t = 0; t < Layout::lattSize()[3]; t++)
          {
            Double a = 0, b=0, b0=0, B=0, B0=0, BO=0, B0O=0, C=0, C0=0, CI=0, C0I=0, c=0, d=0, e=0, e0=0, f=0, f0=0, g=0, g0=0;
            multi1d<int> tCoords;
            tCoords.resize(Nd);
            tCoords[3] = t;
            for(int x = 0; x < Layout::lattSize()[0];x++)
              for(int y = 0; y < Layout::lattSize()[1];y++)
                for(int Z = 0; Z < Layout::lattSize()[2];Z++)
                {
                  tCoords[0] = x; tCoords[1] = y; tCoords[2] = Z;
//                  a += real(trace(peekSite(Oa0, tCoords)));
                  b += real(trace(peekSite(O0, tCoords)));
                  b0 += real(trace(peekSite(O00, tCoords)));
                  B += real(trace(peekSite(O1, tCoords)));
                  B0 += real(trace(peekSite(O10, tCoords)));
                  BO += real(trace(peekSite(O2, tCoords)));
                  B0O += real(trace(peekSite(O20, tCoords)));
                  C += real(trace(peekSite(O3, tCoords)));
                  C0 += real(trace(peekSite(O30, tCoords)));
//                  CI += imag(trace(peekSite(OC, tCoords)));
//                  C0I += imag(trace(peekSite(OC0, tCoords)));
                  c += real(trace(peekSite(Oa, tCoords)));
                  d += real(trace(peekSite(Oa0, tCoords)));
                  e += real(trace(peekSite(Od, tCoords)));
                  e0 += real(trace(peekSite(Od0, tCoords)));
                  f += real(trace(peekSite(OC, tCoords)));
                  f0 += real(trace(peekSite(OC0, tCoords)));
                  g += real(trace(peekSite(OD, tCoords)));
                  g0 += real(trace(peekSite(OD0, tCoords)));
                }

            VecO0.push_back(b);
            VecO00.push_back(b0);
            VecO1.push_back(B);
            VecO10.push_back(B0);
            VecO2.push_back(BO);
            VecO20.push_back(B0O);
            VecO3.push_back(C);
            VecO30.push_back(C0);
            VecOa.push_back(c);
            VecOa0.push_back(d);
            VecOd.push_back(e);
            VecOd0.push_back(e0);
            VecOC.push_back(f);
            VecOC0.push_back(f0);
            VecOD.push_back(g);
            VecOD0.push_back(g0);

          }
/*
        std::vector<Double> VecE0, VecE1, VecE2, VecB0, VecB1, VecB2;
        QDPIO::cout << "Finding En,Bn" << std::endl;
        LatticeColorMatrix E0,E1,E2,B0,B1,B2;
        E0 = E[0]*un*En[0]*adj(un);
        E1 = E[1]*un*En[1]*adj(un);
        E2 = E[2]*un*En[2]*adj(un);
        B0 = B[0]*un*Bn[0]*adj(un);
        B1 = B[1]*un*Bn[1]*adj(un);
        B2 = B[2]*un*Bn[2]*adj(un);

        for(int t = 0; t < Layout::lattSize()[3]; t++)
          {
            Double e0 = 0, e1=0 , e2=0 , b0=0 , b1=0 , b2=0;
            multi1d<int> tCoords;
            tCoords.resize(Nd);
            tCoords[3] = t;
            for(int x = 0; x < Layout::lattSize()[0];x++)
              for(int y = 0; y < Layout::lattSize()[1];y++)
                for(int Z = 0; Z < Layout::lattSize()[2];Z++)
                {
                  tCoords[0] = x; tCoords[1] = y; tCoords[2] = Z;
                  e0 += real(trace(peekSite(E0, tCoords)));
                  e1 += real(trace(peekSite(E1, tCoords)));
                  e2 += real(trace(peekSite(E2, tCoords)));
                  b0 += real(trace(peekSite(B0, tCoords)));
                  b1 += real(trace(peekSite(B1, tCoords)));
                  b2 += real(trace(peekSite(B2, tCoords)));

                }

            VecE0.push_back(e0);
            VecE1.push_back(e1);
            VecE2.push_back(e2);
            VecB0.push_back(b0);
            VecB1.push_back(b1);
            VecB2.push_back(b2);
          }

*/
        write(xml_out, "Oa", VecOa);
        for(int t = 0; t < Layout::lattSize()[3]; t++)
          {
       QDPIO::cout <<"OB   "<< z << "  " << n <<"  "<< t <<"  "<< VecO0.at(t) <<"  "<< VecO00.at(t) <<"  "<< VecO1.at(t) <<"  "<< VecO10.at(t) <<"  "<< VecO2.at(t)<<"  "<< VecO20.at(t) <<"  "<< VecO3.at(t)<<"  "<< VecO30.at(t) <<"  "<< VecOa.at(t) <<"  "<< VecOa0.at(t) <<"  "<< VecOd.at(t) <<"  "<< VecOd0.at(t) << std::endl;

          }


        for(int t = 0; t < Layout::lattSize()[3]; t++)
          {
       QDPIO::cout <<"OPP   "<< z << "  " << n <<"  "<< t <<"  "<< VecOC.at(t) <<"  "<< VecOC0.at(t) <<"  "<< VecOD.at(t) <<"  "<< VecOD0.at(t) << std::endl;

          }

/*
        for(int t = 0; t < Layout::lattSize()[3]; t++)
          {
       QDPIO::cout <<"Ob   "<< z << "  " << n <<"  "<< t <<"  "<< VecOb.at(t) <<"  "<< VecOb0.at(t)<<"  "<< VecObO.at(t) <<"  "<< VecOb0O.at(t)<<"  "<< VecOC.at(t) <<"  "<< VecOC0.at(t)<<"  "<< VecOCI.at(t) <<"  "<< VecOC0I.at(t)<<std::endl;

          }
	
*/       
//        QDPIO::cout << "Finding un" << std::endl;

/*
        LatticeColorMatrix u0, un;
        u0 = u[z];
        un = u[z];
        for(int i=1; i < (nn-1); i++)
        {
            u0 =  shift(u0, FORWARD, z);
            un = un*u0;
        }
*/
/*
        LatticeColorMatrix u0, un, u1, nu;
        multi1d<LatticeColorMatrix> uni;
        uni.resize(2*nn-1);
        u0 = u[z];
        un = u[z];
        u1 = shift(u[z],BACKWARD,z);
        nu = adj(u1);
        uni[1]= u[z];
        uni[1+nn-1]= adj(u1);
        for(int i = 2; i < nn; i++)
        {
            u0 =  shift(u0, FORWARD, z);
            u1 =  shift(u1, BACKWARD, z);
            un = un*u0;
            nu = nu*adj(u1);
            uni[i]=un;
            uni[i+nn-1]=nu;
        }
*/
/*
        QDPIO::cout << "Finding OA" << std::endl;
        
 
        for(int i = 0; i < 4; i++)
          {
                   OA += F[i][z]*un*Fn[i][z]*adj(un)/3;

          }
        
       
        std::vector<Double> VecOA;
        QDPIO::cout << "Finding Tr(OA)" << std::endl;
        for(int t = 0; t < Layout::lattSize()[3]; t++)
          {
            Double a = 0, b=0;
            multi1d<int> tCoords;
            tCoords.resize(Nd);
            tCoords[3] = t;
            for(int x = 0; x < Layout::lattSize()[0];x++)
              for(int y = 0; y < Layout::lattSize()[1];y++)
                for(int z = 0; z < Layout::lattSize()[2];z++)
                {
                  tCoords[0] = x; tCoords[1] = y; tCoords[2] = z;
                  a += real(trace(peekSite(OA, tCoords)));
                }

            VecOA.push_back(a);
//            printf("OB_OUTPUT %4d%13.5f\n", t, a);
          }
        write(xml_out, "OA", VecOA);
        for(int t = 0; t < Layout::lattSize()[3]; t++)
          {
               QDPIO::cout << "OB_OUTPUT" <<" "<< t <<" "<< VecOA.at(t) << std::endl;
          }
*/


/*
        QDPIO::cout << "Finding E/B" << std::endl;
	
       
	multi1d<LatticeColorMatrix> E,B,En,Bn;
        multi2d<LatticeColorMatrix> Eni,Bni;

	E.resize(3);
	B.resize(3);
        En.resize(3);
        Bn.resize(3);
        Eni.resize(3,2*nn-1);
        Bni.resize(3,2*nn-1);


        for(int i = 0; i < 3; i++)
          {
            E[i] = F[3][i];
            En[i] = Fn[3][i];
            int j=(i+1)%3;
            int k=(i+2)%3;
            B[i] += leviCh(i,j,k)*F[j][k];
            Bn[i] += leviCh(i,j,k)*Fn[j][k];
          }

        std::vector<Double> VecOa, VecOa0;
        QDPIO::cout << "Finding Oa" << std::endl;
        LatticeColorMatrix Oa, Oa0;
        for(int i = 0; i < 3;i++)
        {   Oa0 += E[i]*E[i]- B[i]*B[i];
            Oa += E[i]*un*En[i]*adj(un)-B[i]*un*Bn[i]*adj(un);
        }
        QDPIO::cout << "Finding Tr(Oa)" << std::endl;
        for(int t = 0; t < Layout::lattSize()[3]; t++)
          {
            Double a = 0, b=0;
            multi1d<int> tCoords;
            tCoords.resize(Nd);
            tCoords[3] = t;
            for(int x = 0; x < Layout::lattSize()[0];x++)
              for(int y = 0; y < Layout::lattSize()[1];y++)
                for(int Z = 0; Z < Layout::lattSize()[2];Z++)
                {
                  tCoords[0] = x; tCoords[1] = y; tCoords[2] = Z;
                  a += real(trace(peekSite(Oa0, tCoords)));
                  b += real(trace(peekSite(Oa, tCoords)));
                }

            VecOa.push_back(b);
            VecOa0.push_back(a);
          }

        write(xml_out, "Oa", VecOa);
        for(int t = 0; t < Layout::lattSize()[3]; t++)
          {
               QDPIO::cout <<"OB   "<< z  <<"  " << VecOa.at(t) <<"  "<< VecOa0.at(t)  << std::endl;
          }
*/
/*                                                        
        for(int ii=1; ii<(nn); ii++)
        {

	for(int i = 0; i < 3; i++)
	  {
            Eni[i][ii]=Fni[3][i][ii];
	    int j=(i+1)%3;
	    int k=(i+2)%3;
            Bni[i][ii] += leviCh(i,j,k)*Fni[j][k][ii]; 
          }

	std::vector<Double> VecOa, VecOa0;
	QDPIO::cout << "Finding Oa" << std::endl;
	LatticeColorMatrix Oa, Oa0;
        Oa0=0;
        Oa=0;
	for(int i = 0; i < 3;i++)
        {   Oa0 += E[i]*E[i]- B[i]*B[i];
            Oa += E[i]*uni[ii]*Eni[i][ii]*adj(uni[ii])-B[i]*uni[ii]*Bni[i][ii]*adj(uni[ii]);
	}	      
	QDPIO::cout << "Finding Tr(Oa)" << std::endl;
	for(int t = 0; t < Layout::lattSize()[3]; t++)
	  {
	    Double a = 0, b=0;
	    multi1d<int> tCoords;
	    tCoords.resize(Nd);
	    tCoords[3] = t;
	    for(int x = 0; x < Layout::lattSize()[0];x++)
	      for(int y = 0; y < Layout::lattSize()[1];y++)
		for(int Z = 0; Z < Layout::lattSize()[2];Z++)
		{
		  tCoords[0] = x; tCoords[1] = y; tCoords[2] = Z;
		  a += real(trace(peekSite(Oa0, tCoords)));
		  b += real(trace(peekSite(Oa, tCoords)));
                }

	    VecOa.push_back(b);
            VecOa0.push_back(a);
	  }

	write(xml_out, "Oa", VecOa);
        for(int t = 0; t < Layout::lattSize()[3]; t++)
          {
               QDPIO::cout <<"OB   "<< z <<"nn "<< ii <<"  " << VecOa.at(t) <<"  "<< VecOa0.at(t)  << std::endl;
          }
	
       
        }
*/

	// Calculating O_b
	std::vector<Double> O_b;
//	std::vector<Double> E2;
//	std::vector<Double> B2;

        for(int i = 0; i < params.srcs.size(); i++)
	{
	    //Start at t_start
	    QDPIO::cout << "Processing src " << i << " from "
			<< params.srcs[i].t_start << " to "
			<< params.srcs[i].t_end << std::endl;
	    getO_b(O_b, params.srcs[i], params.radius, Oa0);
	}
	    
	//Write O_b out to the xml file
	write(xml_out, "O_b", O_b);
//	write(xml_out, "E2", E2);
//        write(xml_out, "B2", B2);

/*        for(int t = 0; t < Layout::lattSize()[3]; t++)
          {
                      QDPIO::cout <<"OB   "<< z << "  " << n <<"  "<< t <<"  "<< O_b.at(t)<<std::endl;

          }
*/
      }
      }
	pop(xml_out);
		
	snoop.stop();
	QDPIO::cout << InlineObEnv::name << ": total time = "
		    << snoop.getTimeInSeconds() 
		    << " secs" << std::endl;
	
	QDPIO::cout << InlineObEnv::name << ": ran successfully" << std::endl;
	
	END_CODE();

      
      
    } //End func()
/*
    ColorMatrix InlineMyMeas::get_G(const multi1d<int>& coords, int mu, int nu,
				    const multi2d<LatticeColorMatrix>& P)
    {
	//TODO impliment boundry conditions
	ColorMatrix G;
	//Should pull in the four plaquettes
	for(int i = 0; i <= 1; i++)
	    for(int j = 0; j <= 1; j++)
	    {
		multi1d<int> tCoords = coords;
		tCoords[mu] -= i;
		tCoords[nu] -= j;
		G += peekSite(P[mu][nu], tCoords);
		QDPIO::cout << "Getting P_{" << mu << nu << "} at ( ";
		for(int x = 0; x < Nd; x++)
		    QDPIO::cout << tCoords[x] << " ";
		QDPIO::cout << ")" << std::endl;
	    }
	
	return G;
    }
*/





    void InlineMyMeas::getO_b(std::vector<Double>& vecOb,
                              const  InlineObEnv::InlineObParams::Src_t src,
                              const int radius,
                              const LatticeColorMatrix& Oa0)
    {
        multi1d<int> t_coords; //Coords to find O_b at
        t_coords.resize(Nd);
        Double Beta = 1;
        Double a = 1;
        for(int t = src.t_start; t != (src.t_end+1)% Layout::lattSize()[3]; t = (t+1)% Layout::lattSize()[3] )
        {
            Double O_b = 0;
            t_coords[3] = t;
            for(int x = 0; x < Layout::lattSize()[0]; x++)
            {
                t_coords[0] = x;
                for(int y = 0; y < Layout::lattSize()[1]; y++)
                {
                    t_coords[1] = y;
                    for(int z = 0; z < Layout::lattSize()[2]; z++)
                    {
                        t_coords[2] = z;
                        if(validLocation(t_coords, src.srcLoc, radius))
                        {
                            Double e,b;
                            O_b += real(trace(peekSite(Oa0, t_coords)));
                        }
                    }
                }
            } //End sum over space
            vecOb.push_back(O_b);
//        QDPIO::cout <<"sum over space   "<< O_b <<std::endl;
        } //end loop through time
    } //end getO_b

  
/*
    void InlineMyMeas::getO_b(std::vector<Double>& vecOb,
			      std::vector<Double>& vecE,
			      std::vector<Double>& vecB,
			      const  InlineObEnv::InlineObParams::Src_t src,
			      const int radius,
			      const multi2d<LatticeColorMatrix>& plane_plaq)
    {
	//Constants for finding O_b
	multi1d<int> t_coords; //Coords to find O_b at
	t_coords.resize(Nd);
	Double Beta = 1;
	Double a = 1;
	
	for(int t = src.t_start; t != (src.t_end+1)% Layout::lattSize()[3]; t = (t+1)% Layout::lattSize()[3] )
	{
	    //QDPIO::cout << "Processing t=" << t << std::endl;
	    Double O_b = 0;
	    Double E = 0;
	    Double B = 0;
	    t_coords[3] = t;

	    // Sum over all space 
	    for(int x = 0; x < Layout::lattSize()[0]; x++)
	    {
		t_coords[0] = x;
		for(int y = 0; y < Layout::lattSize()[1]; y++)
		{
		    t_coords[1] = y;
		    for(int z = 0; z < Layout::lattSize()[2]; z++)
		    {
			t_coords[2] = z;
			if(validLocation(t_coords, src.srcLoc, radius))
			{
			    Double e,b;
			    O_b += getO_b(t_coords, plane_plaq, e, b);
			    E += e;
			    B += b;
			}
		    } 
		}
	    } //End sum over space

	    // add scale factors to O_b and push back to
	    //   resultant vector 
	    
	    //TODO: Should this be scaled by 1/Nc like the normalized plaquettes are?
	    //O_b *= -1 * Double(4)/9.0 * Beta/a * Double(1)/Nc;
	    vecOb.push_back(O_b);
	    vecE.push_back(E);
	    vecB.push_back(B);
	} //end loop through time
    } //end getO_b
    
    //Code to calculate O_b at coords
    Double InlineMyMeas::getO_b(const multi1d<int>& t_coords,
				const multi2d<LatticeColorMatrix>& plane_plaq,
				Double& E, Double& B)
    {
	E = 0;
	B = 0;
	
	for(int mu = 0; mu < Nd; mu++)
	{
	    if(mu != tDir())
	    {
		//Do first half of sum
		E += real(trace(
				peekSite(plane_plaq[mu][tDir()],
					 t_coords)));
		//do second half of sum
		for(int nu = 0; nu < mu; nu++)
		{
		    if(nu != tDir())
		    {
			B += real(trace(
					peekSite(plane_plaq[nu][mu],
						 t_coords)));
		    }
		}//end second half
	    } //end check for tDir
	} //end 1st half
        return E-B;
    } //end getO_b
*/
    bool InlineMyMeas::validLocation(const multi1d<int>& t_coords,
				     const multi1d<int>& t_src,
				     int R)
    {
	if (R == 0) return true;
	
	int dist = 0;
	for(int i = 0; i < 3; i++)
	{
	    int dx = std::abs(t_coords[i] - t_src[i]);
	    int dimSize = Layout::lattSize()[i];
	    dx = (dx > (dimSize - dx)) ? dimSize - dx : dx;
	    dist += dx*dx;
	}
	
	return dist <= R*R;
	    
    }

};
