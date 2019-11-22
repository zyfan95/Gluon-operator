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
    
    /* This function is used to obtain the gluon field F_\mu\nu shift "len" times along "dir" direction.
     * The input F is one matrix element of F_\mu\nu.
     */
    LatticeColorMatrix  field(int dir, int len, LatticeColorMatrix F)
    {
        LatticeColorMatrix Ftem, Fshift; /* Ftem is the temporary variable,
                                            Fshift is the gluon field F_\mu\nu shift "len" times along "dir" direction that will be returned.
                                          */
        Fshift = F;
        for(int i = 1; i <= len; i++)
        {
                Ftem = Fshift;
                Fshift = shift(Ftem, FORWARD, dir);
        }
        return Fshift;
    }
 
    /* This function is used to construct and return the wilson line at "dir" direction with the length "len".
     * The input u is the unit wilson link at "dir" direction.
     * U(x,x+len\hat{dir})= \sum_{y=x}^{x+len-1}U(y,y+1\hat{dir})
     */
    LatticeColorMatrix  wilsonline(int dir, int len, LatticeColorMatrix u) 
    {
	LatticeColorMatrix utem, ushift, uline; /* utem is the temporary wilson line
						ushift is the shift "len" times wilson link
						uline is the wilson line at "dir" direction with the length "len" that will be returned
						*/
	ushift = u;
	uline = u;
	for(int i = 1; i < len; i++)
        {
		utem = ushift;
		ushift = shift(utem, FORWARD, dir);
		utem = uline;
		uline = utem * ushift;
	}
	return uline; 
    }

    /* The function is to calculte the m*n plaquette at the direction of
     * \mu and \nu. The unit wilson link umu and unu at \mu and \nu direction are
     * necessary input.
     * P1_{\mu\nu}=U(x,x+m\hat{\mu})U(x+m\hat{\mu},x+m\hat{\mu}+n\hat{\nu})U(x+m\hat{\mu}+n\hat{\nu},x+n\hat{\nu})U(x+n\hat{\nu},x)
     * P2_{\mu\nu}=P1_{\nu,-\mu}
     * P3_{\mu\nu}=P1_{-\mu,-\nu}
     * P4_{\mu\nu}=P1_{-\nu,\mu}
     */ 
    LatticeColorMatrix  plaquette(int mu, int nu, int m, int n, LatticeColorMatrix umu, LatticeColorMatrix unu, int pla_num)
    {
        LatticeColorMatrix wltem, wlmu, wlnu, wlmu_shift, wlnu_shift, plane_plaq_mn; /* wltem is the intermediated line
											wlmu is the first quater of the wilson line of the plaquette
                                                                                        wlnu is the second quater of the wilson line of the plaquette
                                                                                        wlmu_shift is the third quater of the wilson line of the plaquette
                                                                                        wlnu_shift is the last quater of the wilson line of the plaquette
											plane_plaq_mn is the m*n plaquette that will be returned
											*/
	wlmu = wilsonline(mu, m, umu); 
	wlnu = wilsonline(nu, n, unu);
	wlmu_shift = wlmu;
        wlnu_shift = wlnu;
        for(int i = 1; i <= n; i++)
        {
                wltem = wlmu_shift;
                wlmu_shift = shift(wltem, FORWARD, nu);
        }
        for(int j = 1; j <= m; j++)
        {
                wltem = wlnu_shift;
                wlnu_shift = shift(wltem, FORWARD, mu);
        }
	switch (pla_num)  //generate 4 plaquette used in the construction of F_\mu\nu
        {
                case 1:plane_plaq_mn = wlmu*wlnu_shift*adj(wlmu_shift)*adj(wlnu);
                break;
                case 2:plane_plaq_mn = shift(wlnu_shift*adj(wlmu_shift)*adj(wlnu)*wlmu, BACKWARD, mu);
                break;
                case 3:plane_plaq_mn = shift(shift(adj(wlmu_shift)*adj(wlnu)*wlmu*wlnu_shift, BACKWARD, mu), BACKWARD, nu);
                break;
                case 4:plane_plaq_mn = shift(adj(wlnu)*wlmu*wlnu_shift*adj(wlmu_shift), BACKWARD, nu);
                break;
        }
	//plane_plaq_mn = ulinem*ulinenshift*adj(ulinemshift)*adj(ulinen);
	return plane_plaq_mn;
    }

    /* The function is used to calculate the mean(real(trace())) of the plaquette.
     * The input should be a matrix element of plane plaquette.
     * The output is a matrix element of the mean(real(trace())) of the plaquette.
     */ 
    Double tr_plane_pla(LatticeColorMatrix plane_plaq)
    {
	Double tr_plane_plaq;
	tr_plane_plaq = sum(real(trace(plane_plaq)));
	tr_plane_plaq /= Double(Layout::vol() * Nc); //average on sites and colors
	tr_plane_plaq = tr_plane_plaq; //symmetric
	return tr_plane_plaq;
    }

    /* The function is used to construct the operators (both local and non-local operators). 
     * The inputs are the wilson line length "len" and the direction "dir", gluon field F_\mu\nu, and the unit wilson link u
     * The outputs are the 6 operators.
     */ 
    multi1d<LatticeColorMatrix> fun_Operator(int len, int dir, multi2d<LatticeColorMatrix> F, LatticeColorMatrix u)
    {
	LatticeColorMatrix un;
	multi2d<LatticeColorMatrix> Fn;
	Fn.resize(Nd,Nd);
	Fn = 0;
	multi1d<LatticeColorMatrix> Op;
	Op.resize(6);
	Op = 0;

	un = wilsonline(dir, len, u); //calculate the wilson line

	for(int mu = 0; mu < Nd; mu++)
        {
        	for(int nu = mu+1; nu < Nd; nu++)
                {
                	Fn[nu][mu] = field(dir, len, F[nu][mu]); //shift the F_\mu\nu
                        Fn[mu][nu] = -Fn[nu][mu];  //anti-symmetric
                }
        }

/*
        LatticeColorMatrix um;
	un = adj(un);
        multi2d<LatticeColorMatrix> Fm;
        Fm.resize(Nd,Nd);
	for(int k=0;k<n; k++)
	{
		um = un;
		un = shift(um, BACKWARD, z);
	}
	for(int mu = 0; mu < Nd; mu++)
        {
                for(int nu = mu+1; nu < Nd; nu++)
                {
			Fm[nu][mu] = F[nu][mu];
                        for(int k=0;k<n; k++)
			{
				Fn[nu][mu] = shift(Fm[nu][mu], BACKWARD, z);
				Fm[nu][mu] = Fn[nu][mu];
			}
                        Fn[mu][nu] = -Fn[nu][mu];  //anti-symmetric
                }
        }
*/	

	/* Operators definition
 	 * O(F_{\mu \nu},F^{\alpha \beta},len) = tr(F_{\mu \nu}(0)U(0,len)F_{\mu \nu}(len)U(0,len))
 	 * O_0 = O(F_{\mu t},F^{\mu t},len)-1/4 O(F_{\mu \nu},F^{\mu \mu},len)
 	 * O_1 = O(F_{\mu t},F^{\mu z},len)
 	 * O_2 = O(F_{\mu z},F^{\mu z},len)-1/4 O(F_{\mu \nu},F^{\mu \mu},len)
 	 * O_3 = O(F_{\mu z},F^{\mu z},len)
 	 * O_4 = O(F_{\mu t},F^{\mu t},len)
 	 * O_d = O(F_{t z},F^{t z},len)
 	 */ 	
	if(len==0) //local operators
	{
		for(int i=0;i<4; i++)
                {
			Op[0] += F[3][i]*F[3][i];
			Op[1] += F[3][i]*F[dir][i];
			Op[2] += F[dir][i]*F[dir][i];
                	Op[3] += F[dir][i]*F[dir][i];
			Op[4] += F[3][i]*F[3][i];
			if(i==3)
			{ Op[5] += F[dir][i]*F[dir][i];}
                }
		for(int i=0;i<4; i++)
                {
                        for(int j=0;j<i;j++)
                        {
                                Op[0] -= 0.5*F[j][i]*F[j][i];
                                Op[2] -= 0.5*F[j][i]*F[j][i];
                        }
               	}

 
	}
	else //non-local
	{
		for(int i=0;i<4; i++)
                {
                        Op[0] += F[3][i]*un*Fn[3][i]*adj(un);
                        Op[1] += F[3][i]*un*Fn[dir][i]*adj(un);
                        Op[2] += F[dir][i]*un*Fn[dir][i]*adj(un);
                        Op[3] += F[dir][i]*un*Fn[dir][i]*adj(un);
                        Op[4] += F[3][i]*un*Fn[3][i]*adj(un);
                        if(i==dir || i==3)
                        { Op[5] += F[dir][i]*un*Fn[dir][i]*adj(un);}
                }
                for(int i=0;i<4; i++)
                {
                        for(int j=0;j<i;j++)
                        {
                                Op[0] -= 0.5*F[j][i]*un*Fn[j][i]*adj(un);
                                Op[2] -= 0.5*F[j][i]*un*Fn[j][i]*adj(un);
                        }
                }
/*
                for(int i=0;i<4; i++)
                {
                        Op[0] += Fn[3][i]*adj(un)*F[3][i]*un;
                        Op[1] += Fn[3][i]*adj(un)*F[z][i]*un;
                        Op[2] += Fn[z][i]*adj(un)*F[z][i]*un;
                        Op[3] += Fn[z][i]*adj(un)*F[z][i]*un;
                        Op[4] += Fn[3][i]*adj(un)*F[3][i]*un;
                        if(i==z || i==3)
                        { Op[5] += Fn[z][i]*adj(un)*F[z][i]*un;}
                }
                for(int i=0;i<4; i++)
                {
                        for(int j=0;j<i;j++)
                        {
                                Op[0] -= 0.5*Fn[j][i]*adj(un)*F[j][i]*un;
                                Op[2] -= 0.5*Fn[j][i]*adj(un)*F[j][i]*un;
                        }
                }
*/


	}

	return Op;
	
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
	
	//chroma function to calculate the plaquette
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
		plane_plaq[0][mu][nu] = plane_plaq[0][nu][mu];
	    }
	}
        
	multi2d<LatticeColorMatrix> plane_plaq_11, plane_plaq_12, plane_plaq_13, plane_plaq_14, plane_plaq_22, plane_plaq_32; //plane_plaq_{m}{n} represents the n*m plaquette
        multi2d<Double> tr_plane_plaq_11, tr_plane_plaq_12, tr_plane_plaq_13, tr_plane_plaq_14, tr_plane_plaq_22, tr_plane_plaq_32; //The mean(trace(real of the plaquette
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
                }
        }
	for(int mu = 0; mu < Nd; mu++)
	{
		for(int nu = mu+1; nu < Nd; nu++)
		{
			tr_plane_plaq_11[mu][nu] = tr_plane_pla(plane_plaq_11[nu][mu]);
                        tr_plane_plaq_12[mu][nu] = tr_plane_pla(plane_plaq_12[nu][mu]);
                        tr_plane_plaq_13[mu][nu] = tr_plane_pla(plane_plaq_13[nu][mu]);
                        tr_plane_plaq_14[mu][nu] = tr_plane_pla(plane_plaq_14[nu][mu]);
                        tr_plane_plaq_22[mu][nu] = tr_plane_pla(plane_plaq_22[nu][mu]);
                        tr_plane_plaq_32[mu][nu] = tr_plane_pla(plane_plaq_32[nu][mu]);
		}
	}	
	
	//Pirnt the mean(trace(real of the plaquette.
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

	//loop over direction 0,1,2
	for(int dir = 0; dir < 3; dir++)
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
		
		// Obtain the F_21 at t=8, y=12, z=12 along x direction
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

		// Print the F_\mu\nu at t=8, y=12, z=12 along x direction
                for(int x = 0; x < Layout::lattSize()[0]; x++)
                {
                        QDPIO::cout <<"F21   "<< x << "   " << VecFmnre.at(x) << "   " << VecFmnim.at(x) << std::endl;

                }

		// Obtain and print the sum(trace( of F_\mu\nu to check the anti-symmetric.
                for(int t = 0; t < Layout::lattSize()[3]; t++)
                {
                	multi2d<Double> F_sumre, F_sumim; //real part and imaginary part of the sum(trace(F_\mu\nu))
			F_sumre.resize(Nd,Nd);
			F_sumim.resize(Nd,Nd);
			F_sumre = 0;
                        F_sumim = 0;
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
				tCoords[2] = Z;
				for(int mu = 0; mu < Nd; mu++)
                		{
                        	 for(int nu = 0; nu < Nd; nu++)
                        	 {
					F_sumre[mu][nu] += real(trace(peekSite(F[mu][nu], tCoords)));
					F_sumim[mu][nu] += imag(trace(peekSite(F[mu][nu], tCoords)));
				 }
				}	
			  }	
			 }
			}
                        for(int mu = 0; mu < Nd; mu++)
                        {
                        	for(int nu = 0; nu < Nd; nu++)
                        	{
					QDPIO::cout <<"F   "<< t << "   "<< mu << "   " << nu << "   "  << F_sumre[mu][nu] << "   " << F_sumim[mu][nu] << std::endl;
				}
			}
          	}
		// loop over the wilson line length len
		for(int len = 0; len < mn; len++)
        	{
			multi1d<LatticeColorMatrix> Op;
			Op.resize(6);
			Op = 0;
			Op = fun_Operator(len, dir, F, u[dir]); //Obtain the operators
			//Calculate and print out the real(trace of the operators
			for(int t = 0; t < Layout::lattSize()[3]; t++)
                        {
                                multi1d<Double> op;
				op.resize(6);
				op = 0;
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
						for(int i = 0; i < 6;i++)
						{
                                                	op[i] += real(trace(peekSite(Op[i], tCoords)));
						}
					}
				QDPIO::cout <<"Op   "<< dir << "  " << len <<"  "<< t <<"  "<< op[0] <<"  "<< op[1] <<"  "<< op[2] <<"  "<< op[3] <<"  "<< op[4]<<"  "<< op[5]<< std::endl;
				
			}//end of t loop
		}//end of len loop
        



	}//end of dir loop


    } 

}
                   
