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


    LatticeColorMatrix  plaquette(int mu, int nu, int m, int n, LatticeColorMatrix umu, LatticeColorMatrix unu, int pla_num) /* The function is to calculte the m*n plaquette at the direction of
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
	switch (pla_num)  //generate 4 plaquette used in the construction of F_\mu\nu
        {
                case 1:plane_plaq_mn = ulinem*ulinenshift*adj(ulinemshift)*adj(ulinen);
                break;
                case 2:plane_plaq_mn = shift(ulinenshift*adj(ulinemshift)*adj(ulinen)*ulinem, BACKWARD, mu);
                break;
                case 3:plane_plaq_mn = shift(shift(adj(ulinemshift)*adj(ulinen)*ulinem*ulinenshift, BACKWARD, mu), BACKWARD, nu);
                break;
                case 4:plane_plaq_mn = shift(adj(ulinen)*ulinem*ulinenshift*adj(ulinemshift), BACKWARD, nu);
                break;
        }
	//plane_plaq_mn = ulinem*ulinenshift*adj(ulinemshift)*adj(ulinen);
	return plane_plaq_mn;
    }

    Double tr_plane_pla(int m, int n, LatticeColorMatrix plane_plaq)
    {
	Double tr_plane_plaq;
	tr_plane_plaq = sum(real(trace(plane_plaq)));
	tr_plane_plaq /= Double(Layout::vol() * Nc);
	tr_plane_plaq = tr_plane_plaq; //symmetric
	return tr_plane_plaq;
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
        plane_plaq.resize(4,Nd,Nd); //Quadrent, plane

        /* Calculate the plaquette and the average trace of
         * the plaquette
         */
        for(int mu = 0; mu < Nd; mu++)
        {
            for(int nu = mu+1; nu < Nd; nu++)
            {
                plane_plaq[0][nu][mu] = plaquette(mu, nu, 1, 1, u[mu], u[nu], 1);
		plane_plaq[1][nu][mu] = plaquette(mu, nu, 1, 1, u[mu], u[nu], 2);
                plane_plaq[2][nu][mu] = plaquette(mu, nu, 1, 1, u[mu], u[nu], 3);
                plane_plaq[3][nu][mu] = plaquette(mu, nu, 1, 1, u[mu], u[nu], 4);
	    }
	}
        
	multi2d<LatticeColorMatrix> plane_plaq_11, plane_plaq_12, plane_plaq_13, plane_plaq_14, plane_plaq_22, plane_plaq_32;
        multi2d<Double> tr_plane_plaq_11, tr_plane_plaq_12, tr_plane_plaq_13, tr_plane_plaq_14, tr_plane_plaq_22, tr_plane_plaq_32;
        plane_plaq_11.resize(Nd,Nd);
        plane_plaq_12.resize(Nd,Nd);
        plane_plaq_13.resize(Nd,Nd);
        plane_plaq_14.resize(Nd,Nd);
        plane_plaq_22.resize(Nd,Nd);
        plane_plaq_32.resize(Nd,Nd);
        tr_plane_plaq_11.resize(Nd,Nd);
        tr_plane_plaq_12.resize(Nd,Nd);
        tr_plane_plaq_13.resize(Nd,Nd);
        tr_plane_plaq_14.resize(Nd,Nd);
        tr_plane_plaq_22.resize(Nd,Nd);
        tr_plane_plaq_32.resize(Nd,Nd);

        for(int mu = 0; mu < Nd; mu++)
        {
        	for(int nu = mu+1; nu < Nd; nu++)
                {
			plane_plaq_11[nu][mu] = plaquette(mu, nu, 1, 1, u[mu], u[nu], 1);
			plane_plaq_12[nu][mu] = plaquette(mu, nu, 1, 2, u[mu], u[nu], 1);
                       	plane_plaq_13[nu][mu] = plaquette(mu, nu, 1, 3, u[mu], u[nu], 1);
                        plane_plaq_14[nu][mu] = plaquette(mu, nu, 1, 4, u[mu], u[nu], 1);
                       	plane_plaq_22[nu][mu] = plaquette(mu, nu, 2, 2, u[mu], u[nu], 1);
                        plane_plaq_32[nu][mu] = plaquette(mu, nu, 3, 2, u[mu], u[nu], 1);
			plane_plaq_11[nu][mu] = plane_plaq_11[mu][nu];
			plane_plaq_12[nu][mu] = plane_plaq_12[mu][nu];
                        plane_plaq_13[nu][mu] = plane_plaq_13[mu][nu];
                        plane_plaq_14[nu][mu] = plane_plaq_14[mu][nu];
                        plane_plaq_22[nu][mu] = plane_plaq_22[mu][nu];
                        plane_plaq_32[nu][mu] = plane_plaq_32[mu][nu];
                }
        }
	for(int mu = 0; mu < Nd; mu++)
	{
		for(int nu = mu+1; nu < Nd; nu++)
		{
			tr_plane_plaq_11[mu][nu] = tr_plane_pla(mu, nu, plane_plaq_11[nu][mu]); //symmetric
                        tr_plane_plaq_12[mu][nu] = tr_plane_pla(mu, nu, plane_plaq_12[nu][mu]); //symmetric
                        tr_plane_plaq_13[mu][nu] = tr_plane_pla(mu, nu, plane_plaq_13[nu][mu]); //symmetric
                        tr_plane_plaq_14[mu][nu] = tr_plane_pla(mu, nu, plane_plaq_14[nu][mu]); //symmetric
                        tr_plane_plaq_22[mu][nu] = tr_plane_pla(mu, nu, plane_plaq_22[nu][mu]); //symmetric
                        tr_plane_plaq_32[mu][nu] = tr_plane_pla(mu, nu, plane_plaq_32[nu][mu]); //symmetric
		}
	}	

        for(int mu = 0; mu < Nd; mu++)
        {
        	for(int nu = mu+1; nu < Nd; nu++)
        	{
                	write(xml_out, "plane_plaq_11_" + std::to_string(mu) +
                	std::to_string(nu), tr_plane_plaq_11[mu][nu]);
                }
        }

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
                        write(xml_out, "plane_plaq_32_" + std::to_string(mu) +
                        std::to_string(nu), tr_plane_plaq_32[mu][nu]);
                }
        }


	for(int z = 0; z < 3; z++)
        {
                QDPIO::cout << "Finding F" << std::endl;
                int mn = 6; //maximum wilson length
                /** Find F_{n,munu} **/
                multi2d<LatticeColorMatrix> F, Fn;
                F.resize(Nd,Nd);
                Fn.resize(Nd,Nd);
		F = 0;
		Fn = 0;

                // Calculate the F
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
                        }
                 }
                std::vector<Double> VecFmnre, VecFmnim;
                for(int x = 0; x < Layout::lattSize()[0]; x++)
                {
                        Double fmnre = 0, fmnim = 0;
                        multi1d<int> fCoords;
                        fCoords.resize(Nd);
                        fCoords[3] = 8;
                        fCoords[0] = x;
                        fCoords[1] = 12;
                        fCoords[2] = 12;
                        fmnre = real(trace(peekSite(F[2][1], fCoords)));
                        fmnim = imag(trace(peekSite(F[2][1], fCoords)));
                        VecFmnre.push_back(fmnre);
			VecFmnim.push_back(fmnim);
                }

                for(int x = 0; x < Layout::lattSize()[0]; x++)
                {
                        QDPIO::cout <<"reF21   "<< x << "   " << VecFmnre.at(x) << "   " << VecFmnim.at(x) << std::endl;

                }

                std::vector<Double> VecFmeanre, VecFmeanim;
                for(int t = 0; t < Layout::lattSize()[3]; t++)
                {
                	multi2d<Double> Fmeanre, Fmeanim;
			Fmeanre.resize(Nd,Nd);
			Fmeanim.resize(Nd,Nd);
			Fmeanre = 0;
                        Fmeanim = 0;
                 	multi1d<int> tCoords;
                        tCoords.resize(Nd);
                        tCoords[3] = t;
                	for(int x = 0; x < Layout::lattSize()[0]; x++)
                	{
			 for(int y = 0; y < Layout::lattSize()[0]; y++)
                	 {
                          for(int Z = 0; Z < Layout::lattSize()[0]; Z++)
                          {
			  	tCoords[0] = x;
				tCoords[1] = y;
				tCoords[2] = z;
				for(int mu = 0; mu < Nd; mu++)
                		{
                        	 for(int nu = 0; nu < Nd; nu++)
                        	 {
					Fmeanre[mu][nu] += real(trace(peekSite(F[mu][nu], tCoords)));
					Fmeanim[mu][nu] += imag(trace(peekSite(F[mu][nu], tCoords)));
				 }
				}	
			  }	
			 }
			}
                        for(int mu = 0; mu < Nd; mu++)
                        {
                        	for(int nu = 0; nu < Nd; nu++)
                        	{
					QDPIO::cout <<"Fmean   "<< t << "   "<< mu << "   " << nu << "   "  << Fmeanre[mu][nu] << "   " << Fmeanim[mu][nu] << std::endl;
				}
			}
          	}

		for(int n = 1; n < mn; n++)
        	{
			LatticeColorMatrix un;
			un = wilsonline(z, n, u[z]);

                        //Calculate Fn for different n
                        for(int mu = 0; mu < Nd; mu++)
                        {
                                for(int nu = mu+1; nu < Nd; nu++)
                                {
                                        Fn[nu][mu] = field(z, n, F[nu][mu]);
                                        Fn[mu][nu] = -Fn[nu][mu];  //anti-symmetric
                                }
                        }


        		LatticeColorMatrix O0, O00;
        		O0 = 0;


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

		}
        



	}


    } 

}
                   
