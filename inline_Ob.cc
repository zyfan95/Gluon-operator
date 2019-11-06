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
    
   
    LatticeColorMatrix  field(int z, int n, LatticeColorMatrix F) /* This function is used to obtain the gluon field F_\mu\nu shift n times along z direction.
                                                                        The input F is one matrix element of F_\mu\nu.
                                                                        */
    {
        LatticeColorMatrix Fmed, Fshift; /* Fmed is the intermediated variable,
                                                Fshift is the gluon field F_\mu\nu shift n times along z direction that will be returned
                                                */
        Fshift = F;
        for(int i = 1; i <= n; i++)
        {
                Fmed = Fshift;
                Fshift = shift(Fmed, FORWARD, z);
        }
        return Fshift;
    }
 

    LatticeColorMatrix  wilsonline(int z, int n, LatticeColorMatrix u) /* This function is used to construct and return the wilson line at z direction with the length n. 
									The input u is the unit wilson link at z direction. 
									*/
    {
	LatticeColorMatrix umed, ushift, uline; /* umed is the intermediated line
						ushift is the shift n times wilson link
						uline is the wilson line at z direction with the length n that will be returned
						*/
	ushift = u;
	uline = u;
	for(int i = 1; i < n; i++)
        {
		umed = ushift;
		ushift = shift(umed, FORWARD, z);
		umed = uline;
		uline = umed * ushift;
	}
	return uline; 
    }


    LatticeColorMatrix  plaquette(int mu, int nu, int m, int n, LatticeColorMatrix umu, LatticeColorMatrix unu) /* The function is to calculte the m*n plaquette at the direction of
														\mu and \nu. The unit wilson link umu and unu at \mu and \nu direction are
														necessary input. 
														*/
    {
        LatticeColorMatrix umed, ulinem, ulinen, ulinemshift, ulinenshift, plane_plaq_mn; /* umed is the intermediated line
											ulinem is the first quater of the wilson line of the plaquette
                                                                                        ulinen is the second quater of the wilson line of the plaquette
                                                                                        ulinemshift is the third quater of the wilson line of the plaquette
                                                                                        ulinenshift is the last quater of the wilson line of the plaquette
											plane_plaq_mn is the m*n plaquette that will be returned
											*/
	ulinem = wilsonline(mu, m, umu); 
	ulinen = wilsonline(nu, n, unu);
	ulinemshift = ulinem;
        ulinenshift = ulinen;
        for(int i = 1; i <= n; i++)
        {
                umed = ulinemshift;
                ulinemshift = shift(umed, FORWARD, nu);
        }
        for(int j = 1; j <= m; j++)
        {
                umed = ulinenshift;
                ulinenshift = shift(umed, FORWARD, mu);
        }
	plane_plaq_mn = ulinem*ulinenshift*adj(ulinemshift)*adj(ulinen);
	return plane_plaq_mn;
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

	Wloop(xml_out, "GMF_O_b", u);
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
	
        
        /* Calculate the different operators for different direction z and 
         * wilson link length n.
        */

        for(int z = 0; z < 3; z++)
        {
		QDPIO::cout << "Finding F" << std::endl;
        	int mn = 6; //maximum wilson length
		/** Find F_{n,munu} **/
		multi2d<LatticeColorMatrix> F, Fn, nF, Fm;
		F.resize(Nd,Nd);
        	Fn.resize(Nd,Nd);
        	Fm.resize(Nd,Nd);
        	nF.resize(Nd,Nd);
                
                // Initialize the link with n=1 along z direction and backward z directon
        	LatticeColorMatrix u0, un, u1, u2, u3;
        	u0 = u[z];
       	 	un = u[z]; //First forward link
        	//u1 = adj(shift(u[z], BACKWARD, z)); 
        	//nu = u1; //First backward link

		// Calculate the F, set the Fn and nF equal to the local one.
        	for(int mu = 0; mu < Nd; mu++)
        	{
            		for(int nu = mu+1; nu < Nd; nu++)
            		{
                	for(int i = 0; i < 4; i++)
                	{
                    		F[nu][mu] += plane_plaq[i][nu][mu] - adj(plane_plaq[i][nu][mu]);
               		}
                	F[nu][mu] *= 1/8.0;
                	F[mu][nu]= -F[nu][mu]; //anti-symmetric
                	Fn[nu][mu] = F[nu][mu];
                	Fn[mu][nu] = -Fn[nu][mu];
                	nF[nu][mu] = F[nu][mu];
                	nF[mu][nu] = -nF[nu][mu];
           		}	
       		 }   

/*
                        for(int i=1;i<4; i++)
                        {
                                O1 += F[3][i]*un*Fn[z][i]*adj(un);
                                O10 += F[3][i]*F[z][i];
                        }
 
                                multi1d<int> Coord_loc, Coord_shift ;
                                Coord_loc.resize(Nd);
				Coord_shift.resize(Nd);
                                Coord_loc[3] = t;
                                Coord_shift[3] = t;
                                for(int x = 0; x < Layout::lattSize()[0];x++)
                                  for(int y = 0; y < Layout::lattSize()[1];y++)
                                    for(int Z = 0; Z < Layout::lattSize()[2];Z++)
                                        {
                                                Coord_loc[0] = x;
                                                Coord_loc[1] = y;
                                                Coord_loc[2] = Z;
						Coord_shift[0] = x;
                                                Coord_shift[1] = y;
                                                Coord_shift[2] = Z;
						Coord_shift[z] = (Coord_shift[z]+1)%Layout::lattSize()[z];
						 =peekSite(F[3][i], Coord_loc)*peekSite(u[z],Coord_loc)*peekSite(Fn[z][i], Coord_shift)*peekSite(adj(u[z]),Coord_loc);	
					}
*/
/*
		
                multi1d<LatticeColorMatrix> u_t, u_m, u_p, u_v;

                u_t.resize(4);
		u_m.resize(4);
		u_p.resize(4);
		u_v.resize(4);


                for(int mu = 0; mu < Nd; mu++)
                {
                        u_t[mu] = u_m[mu];
                        u_m[mu] = shift(u_t[mu], FORWARD, mu);
                        u_t[mu] = u_p[mu];
                        u_p[mu] = u_t[mu]*u_m[mu];
                        u_t[mu] = u_v[mu];
                        u_v[mu] = ;
                }


                for(int mu = 0; mu < Nd; mu++)
                {
                        for(int nu = mu+1; nu < Nd; nu++)
                        {

*/

                //loop over the wilson length n
        	for(int n = 1; n < mn; n++)
        	{
			//Calculate Fn for different n
			for(int mu = 0; mu < Nd; mu++)
			{
	    			for(int nu = mu+1; nu < Nd; nu++)
	    			{

			                Fm[nu][mu] = Fn[nu][mu];
                		//	Fn[nu][mu] = shift(Fm[nu][mu], FORWARD, z); //shift Fn forward one lattice unit along z direction
                			Fn[nu][mu] = field(z, n, F[nu][mu]);	
                			Fn[mu][nu] = -Fn[nu][mu];  //anti-symmetric

                                        Fm[nu][mu] = nF[nu][mu];
                			nF[nu][mu] = shift(Fm[nu][mu], BACKWARD, z);  //shift nF backward one lattice unit along z direction
                			nF[mu][nu] = -nF[nu][mu];
/*
                			Fn[nu][mu] = shift(Fn[nu][mu], FORWARD, z);
                			Fn[mu][nu] = -Fn[nu][mu];
                			nF[nu][mu] = shift(nF[nu][mu], BACKWARD, z);
                			nF[mu][nu] = -nF[nu][mu];
*/

            			}
        		}
			//calculate the wilson link              
        		if(n!=1)
        		{
          			
          			u2 =  shift(u0, FORWARD, z);
                    		u0 = u2;
				u3 = un;
          			un = u3*u2; //add the nth link to the 1~n-1 link
          			//u1 =  shift(u1, BACKWARD, z);
          			//nu = nu*u1;

        		}		
/*
                multi2d<LatticeColorMatrix> plane_plaq_1n;
                multi2d<Double> tr_plane_plaq_1n;
                plane_plaq_1n.resize(Nd,Nd);
                tr_plane_plaq_1n.resize(Nd,Nd);

                for(int mu = 0; mu < Nd; mu++)
                {
                        u_t[mu] = u_m[mu];
                        u_m[mu] = shift(u_t[mu], FORWARD, mu);
                        u_t[mu] = u_p[mu];
                        u_p[mu] = u_t[mu]*u_m[mu];
			for(int nu = mu+1; nu < Nd; nu++)
                        {
                        	u_t[nu] = u_v[nu];
                        	u_v[nu] = shift(u_t[nu], FORWARD, mu);
			}
                }


                for(int mu = 0; mu < Nd; mu++)
                {
                        for(int nu = mu+1; nu < Nd; nu++)
                        {

*/
		un = wilsonline(z, n, u[z]);

        	multi2d<LatticeColorMatrix> plane_plaq_12, plane_plaq_21, plane_plaq_13, plane_plaq_14, plane_plaq_22, plane_plaq_23;
        	multi2d<Double> tr_plane_plaq_12, tr_plane_plaq_21, tr_plane_plaq_13, tr_plane_plaq_14, tr_plane_plaq_22, tr_plane_plaq_23;
        	plane_plaq_12.resize(Nd,Nd);
		plane_plaq_21.resize(Nd,Nd);
                plane_plaq_13.resize(Nd,Nd);
                plane_plaq_14.resize(Nd,Nd);
                plane_plaq_22.resize(Nd,Nd);
                plane_plaq_23.resize(Nd,Nd);
        	tr_plane_plaq_12.resize(Nd,Nd);
		tr_plane_plaq_21.resize(Nd,Nd);
                tr_plane_plaq_13.resize(Nd,Nd);
                tr_plane_plaq_14.resize(Nd,Nd);
                tr_plane_plaq_22.resize(Nd,Nd);
                tr_plane_plaq_23.resize(Nd,Nd);

		for(int mu = 0; mu < Nd; mu++)
        	{
            		for(int nu = mu+1; nu < Nd; nu++)
            		{
			//	plane_plaq_12[nu][mu] = u[mu]*shift(u[mu], FORWARD, mu)*shift(shift(u[nu], FORWARD, mu), FORWARD, mu)*shift(adj(u[mu]*shift(u[mu], FORWARD, mu)), FORWARD, nu)*adj(u[nu]);
			//	plane_plaq_12[nu][mu] = u[mu]*shift(u[mu], FORWARD, mu)*shift(shift(u[nu], FORWARD, mu), FORWARD, mu)*adj(shift(u[mu]*shift(u[mu], FORWARD, mu), FORWARD, nu))*adj(u[nu]);
				plane_plaq_12[nu][mu] = plaquette(mu, nu, 1, 2, u[mu], u[nu]);
				plane_plaq_21[nu][mu] = plaquette(mu, nu, 2, 1, u[mu], u[nu]);
                                plane_plaq_13[nu][mu] = plaquette(mu, nu, 1, 3, u[mu], u[nu]);
                                plane_plaq_14[nu][mu] = plaquette(mu, nu, 1, 4, u[mu], u[nu]);
                                plane_plaq_22[nu][mu] = plaquette(mu, nu, 2, 2, u[mu], u[nu]);
                                plane_plaq_23[nu][mu] = plaquette(mu, nu, 3, 2, u[mu], u[nu]);
			}
		}
		/*
		if(n == 2)
		{
			 plane_plaq_12[z+1][z] = un*shift(shift(u[z+1], FORWARD, z), FORWARD, z)*shift(adj(un), FORWARD, z+1)*adj(u[z+1]); 			
		}
		*/
                for(int mu = 0; mu < Nd; mu++)
                {
                        for(int nu = mu+1; nu < Nd; nu++)
                        {
                                tr_plane_plaq_12[nu][mu] = sum(real(trace(plane_plaq_12[nu][mu])));
                                tr_plane_plaq_12[nu][mu] /= Double(Layout::vol() * Nc);
                                plane_plaq_12[mu][nu] = plane_plaq_12[nu][mu]; //symmetric
                                tr_plane_plaq_12[mu][nu] = tr_plane_plaq_12[nu][mu]; //symmetric
                                tr_plane_plaq_21[nu][mu] = sum(real(trace(plane_plaq_21[nu][mu])));
                                tr_plane_plaq_21[nu][mu] /= Double(Layout::vol() * Nc);
                                plane_plaq_21[mu][nu] = plane_plaq_21[nu][mu]; //symmetric
                                tr_plane_plaq_21[mu][nu] = tr_plane_plaq_21[nu][mu]; //symmetric

                                tr_plane_plaq_13[nu][mu] = sum(real(trace(plane_plaq_13[nu][mu])));
                                tr_plane_plaq_13[nu][mu] /= Double(Layout::vol() * Nc);
                                plane_plaq_13[mu][nu] = plane_plaq_13[nu][mu]; //symmetric
                                tr_plane_plaq_13[mu][nu] = tr_plane_plaq_13[nu][mu]; //symmetric

                                tr_plane_plaq_14[nu][mu] = sum(real(trace(plane_plaq_14[nu][mu])));
                                tr_plane_plaq_14[nu][mu] /= Double(Layout::vol() * Nc);
                                plane_plaq_14[mu][nu] = plane_plaq_14[nu][mu]; //symmetric
                                tr_plane_plaq_14[mu][nu] = tr_plane_plaq_14[nu][mu]; //symmetric

                                tr_plane_plaq_22[nu][mu] = sum(real(trace(plane_plaq_22[nu][mu])));
                                tr_plane_plaq_22[nu][mu] /= Double(Layout::vol() * Nc);
                                plane_plaq_22[mu][nu] = plane_plaq_22[nu][mu]; //symmetric
                                tr_plane_plaq_22[mu][nu] = tr_plane_plaq_22[nu][mu]; //symmetric

                                tr_plane_plaq_23[nu][mu] = sum(real(trace(plane_plaq_23[nu][mu])));
                                tr_plane_plaq_23[nu][mu] /= Double(Layout::vol() * Nc);
                                plane_plaq_23[mu][nu] = plane_plaq_23[nu][mu]; //symmetric
                                tr_plane_plaq_23[mu][nu] = tr_plane_plaq_23[nu][mu]; //symmetric


                        }
                }





                /** Write plane plaq to xml file **/

        	for(int mu = 0; mu < Nd; mu++)
        	{
            		for(int nu = mu+1; nu < Nd; nu++)
            		{
                		write(xml_out, "plane_plaq_12_" + std::to_string(mu) +
                      		std::to_string(nu), tr_plane_plaq_12[mu][nu]);
            		}
        	}

                for(int mu = 0; mu < Nd; mu++)
                {
                        for(int nu = mu+1; nu < Nd; nu++)
                        {
                                write(xml_out, "plane_plaq_21_" + std::to_string(mu) +
                                std::to_string(nu), tr_plane_plaq_21[mu][nu]);
                        }
                }

                for(int mu = 0; mu < Nd; mu++)
                {
                        for(int nu = mu+1; nu < Nd; nu++)
                        {
                                write(xml_out, "plane_plaq_13_" + std::to_string(mu) +
                                std::to_string(nu), tr_plane_plaq_13[mu][nu]);
                        }
                }

                for(int mu = 0; mu < Nd; mu++)
                {
                        for(int nu = mu+1; nu < Nd; nu++)
                        {
                                write(xml_out, "plane_plaq_14_" + std::to_string(mu) +
                                std::to_string(nu), tr_plane_plaq_14[mu][nu]);
                        }
                }

                for(int mu = 0; mu < Nd; mu++)
                {
                        for(int nu = mu+1; nu < Nd; nu++)
                        {
                                write(xml_out, "plane_plaq_22_" + std::to_string(mu) +
                                std::to_string(nu), tr_plane_plaq_22[mu][nu]);
                        }
                }

                for(int mu = 0; mu < Nd; mu++)
                {
                        for(int nu = mu+1; nu < Nd; nu++)
                        {
                                write(xml_out, "plane_plaq_23_" + std::to_string(mu) +
                                std::to_string(nu), tr_plane_plaq_23[mu][nu]);
                        }
                }

                                                	
/*
		        QDPIO::cout << "Finding E/B" << std::endl;


		        multi1d<LatticeColorMatrix> E,B,En,Bn,nE,nB;

		        E.resize(3);
			B.resize(3);
        		En.resize(3);
        		Bn.resize(3);
        		nE.resize(3);
        		nB.resize(3);
*/
                        QDPIO::cout << "Finding O" << std::endl;


        		LatticeColorMatrix O0, O00, O1, O10, O2, O20, O3, O30, Oa, Oa0, Od, Od0, OC, OC0, OD, OD0;
			
			//Initialize the operators *Need to check!*
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


			//Constucting O0, O1, O2, O3, Oa, Od
        
		        for(int i=0;i<4; i++)
       			{
            			O3 += F[z][i]*un*Fn[z][i]*adj(un); //non-local operator
            			O30 += F[z][i]*F[z][i]; //local operator
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

       			//O1 = F[3][0]*un*Fn[z][0]*adj(un)+F[3][1]*un*Fn[z][1]*adj(un)+F[3][2]*un*Fn[z][2]*adj(un);
        		//O10 = F[3][0]*F[z][0]+F[3][1]*F[z][1]+F[3][2]*F[z][2];

		        for(int i=1;i<4; i++)
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
/*
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
*/
/*        		LatticeColorMatrix Oa, Oam, Oa0, Oc, Oc0, Od, Od0;
        		for(int i = 0; i < 3;i++)
        		{   	//O00 += E[i]*E[i]- B[i]*B[i];
            			//O0 += E[i]*un*En[i]*adj(un)-B[i]*un*Bn[i]*adj(un);
            			Oam += E[i]*nu*nE[i]*adj(nu)-B[i]*nu*nB[i]*adj(nu);

        
            			if(i!=z)
            			{
            				OC0 += E[i]*B[i];
            				OC += E[i]*un*Bn[i]*adj(un);
            				OD0 += E[i]*nu*nB[i]*adj(nu);
            				OD += B[i]*un*En[i]*adj(un);

            			}

        		}
*/
			//Finding Tr(O), save them in the VecO
        		QDPIO::cout << "Finding Tr(O)" << std::endl;
			std::vector<Double> VecO0, VecO00, VecO1, VecO10, VecO2, VecO20, VecO3, VecO30, VecOa, VecOa0, VecOd, VecOd0, VecOC, VecOC0, VecOD, VecOD0;
        		for(int t = 0; t < Layout::lattSize()[3]; t++)
          		{
            			Double a = 0, b=0, b0=0, B=0, B0=0, BO=0, B0O=0, C=0, C0=0, CI=0, C0I=0, d=0, d0=0, e=0, e0=0, f=0, f0=0, g=0, g0=0;
           			multi1d<int> tCoords;
            			tCoords.resize(Nd);
            			tCoords[3] = t;
            			for(int x = 0; x < Layout::lattSize()[0];x++)
              			  for(int y = 0; y < Layout::lattSize()[1];y++)
                		    for(int Z = 0; Z < Layout::lattSize()[2];Z++)
                			{
                  				tCoords[0] = x; 
						tCoords[1] = y; 
						tCoords[2] = Z;

                  				b += real(trace(peekSite(O0, tCoords)));
			                  	b0 += real(trace(peekSite(O00, tCoords)));
                  				B += real(trace(peekSite(O1, tCoords)));
                  				B0 += real(trace(peekSite(O10, tCoords)));
                  				BO += real(trace(peekSite(O2, tCoords)));
                  				B0O += real(trace(peekSite(O20, tCoords)));
                  				C += real(trace(peekSite(O3, tCoords)));
                  				C0 += real(trace(peekSite(O30, tCoords)));
                  				d += real(trace(peekSite(Oa, tCoords)));
                  				d0 += real(trace(peekSite(Oa0, tCoords)));
                  				e += real(trace(peekSite(Od, tCoords)));
                  				e0 += real(trace(peekSite(Od0, tCoords)));
                  				//f += real(trace(peekSite(OC, tCoords)));
                 		 		//f0 += real(trace(peekSite(OC0, tCoords)));
                  				//g += real(trace(peekSite(OD, tCoords)));
                  				//g0 += real(trace(peekSite(OD0, tCoords)));
                			}

            				VecO0.push_back(b);
            				VecO00.push_back(b0);
            				VecO1.push_back(B);
            				VecO10.push_back(B0);
            				VecO2.push_back(BO);
            				VecO20.push_back(B0O);
            				VecO3.push_back(C);
            				VecO30.push_back(C0);
            				VecOa.push_back(d);
            				VecOa0.push_back(d0);
            				VecOd.push_back(e);
            				VecOd0.push_back(e0);
            				//VecOC.push_back(f);
            				//VecOC0.push_back(f0);
            				//VecOD.push_back(g);
            				//VecOD0.push_back(g0);

          		}
                        //print the operator at different z, n, t
        		write(xml_out, "Oa", VecOa);
        		for(int t = 0; t < Layout::lattSize()[3]; t++)
          		{
       				QDPIO::cout <<"OB   "<< z << "  " << n <<"  "<< t <<"  "<< VecO0.at(t) <<"  "<< VecO00.at(t) <<"  "<< VecO1.at(t) <<"  "<< VecO10.at(t) <<"  "<< VecO2.at(t)<<"  "<< VecO20.at(t) <<"  "<< VecO3.at(t)<<"  "<< VecO30.at(t) <<"  "<< VecOa.at(t) <<"  "<< VecOa0.at(t) <<"  "<< VecOd.at(t) <<"  "<< VecOd0.at(t) << std::endl;

          		}	

/*
        for(int t = 0; t < Layout::lattSize()[3]; t++)
          {
       QDPIO::cout <<"OPP   "<< z << "  " << n <<"  "<< t <<"  "<< VecOC.at(t) <<"  "<< VecOC0.at(t) <<"  "<< VecOD.at(t) <<"  "<< VecOD0.at(t) << std::endl;

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


//A lot of codes commented below were planning to do cluster decomposition

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
