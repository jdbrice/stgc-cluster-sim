

#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TGraph.h"

#include <iostream>
#include <map>

#define LOGURU_IMPLEMENTATION 1
#include "loguru.h"


const float pent_base = 537.0; // in mm
const float pent_nib = 179.0; // in mm
const float pent_shift = 101.6; // in mm
const int n_clusters_to_gen = 10; // Number of clusters in this "event"
const float noise_prob = 0.2; // 0 -1, 1 is max noise
const float noise_level = 1; // in ADC units (pre integration)
const float strip_pitch = 3.2; // digitize strip pitch
const size_t sat_above = 1024; // ADC
const float cluster_max_adc = 30; // ADC (pre integration)

TF1 * fClusterProfile = nullptr;

using namespace std;
map<string, TH1*> hist;

TGraph * hit_graph = new TGraph();

bool in_bounds( float x, float y ){
    if ( x>0 && x < pent_nib && y < pent_base && y > 0 )
        return true;
    if ( x>0 && x < pent_base && y < pent_nib && y>0 )
        return true;

    float dy = pent_base - pent_nib;
    float dx = pent_nib - pent_base;
    float yedge = (dy/dx) * ( x - pent_nib ) + pent_base;
    // LOG_F( INFO, "yedge = %f @ x = %f", yedge, x );

    if ( y < yedge && y < pent_base && x < pent_base && x>0 && y>0 ) return true;

    return false;
}



int in_bounds_quad( float x, float y ){
    if ( in_bounds( x, y ) ){ // top right
        return 1;
    } else if ( in_bounds( -x, y ) ){ //top left
        return 2;
    } else if ( in_bounds( -x - pent_shift, -y ) ){ // bottom left
        return 3;
    } else if ( in_bounds( x - pent_shift, -y ) ){ // bottom right
        return 4;
    }

    return 0;
}

int map_to_strip( float x, float y, int &sx, int &sy ){
    int q = in_bounds_quad( x, y );
    if ( q == 1 || q == 2 ){
        sx = x / strip_pitch;
        sy = y / strip_pitch;
    } else if ( q == 3 ){
        sx = (x ) / strip_pitch;
        sy = y / strip_pitch;
    } else if ( q == 4 ){
        sx = (x ) / strip_pitch;
        sy = y / strip_pitch;
    }
    return q;
}

int map_to_local_strip( float x, float y, int &sx, int &sy ){
    int q = in_bounds_quad( x, y );
    if ( q == 1 || q == 2 ){
        sx = abs(x / strip_pitch);
        sy = abs(y / strip_pitch);
    } else if ( q == 3 ){
        sx = abs((x + pent_shift ) / strip_pitch);
        sy = abs(y / strip_pitch);
    } else if ( q == 4 ){
        sx = abs((x - pent_shift) / strip_pitch);
        sy = abs(y / strip_pitch);
    }
    return q;
}

int strip_group( int sx, int sy ){
    if ( (sx < 55 && sy < 150) || sy >= 150 ) return 0;
    if ( (sx < 110 && sy < 95) || sy >= 95 ) return 1;
    return 2;
}

void fill_cluster_GEN( float x, float y ){
    TH2 * hGEN = (TH2*)hist[ "hGEN0" ];
    TH2 * hDIG = (TH2*)hist[ "hDIG1" ];
    int bx = hGEN->GetXaxis()->FindBin( x );
    int by = hGEN->GetYaxis()->FindBin( y );

    float maxADC = gRandom->Rndm() * cluster_max_adc;
    LOG_F( INFO, "Cluster at (%f, %f) w maxADC = %f", x, y, maxADC );

    
    int wSpan = 20;
    for ( int i = bx - wSpan; i < bx + wSpan; i++ ){
        for ( int j = by - wSpan; j < by + wSpan; j++ ){
            float bcx = hGEN->GetXaxis()->GetBinCenter( i );
            float bcy = hGEN->GetYaxis()->GetBinCenter( j );

            int q = in_bounds_quad( bcx, bcy );
            if ( q <= 0) continue;

            float r = sqrt( pow( bcx - x, 2 ) + pow( bcy - y, 2 ) );
            float v = fClusterProfile->Eval( r ) * maxADC;
            // LOG_F( INFO, "---- Cluster ADC = %f @ r = %f", v, r );
            hGEN->Fill( bcx, bcy, v );

            int sx = -1, sy = -1;
            map_to_strip( bcx, bcy, sx, sy );
            hDIG->Fill( sx, sy, v );

            map_to_local_strip( bcx, bcy, sx, sy );
            string hname = TString::Format( "hDIG1L%d", q ).Data();
            int xG = strip_group( sx, sy );
            int yG = strip_group( sy, sx ); // not a typo, just reuse function with rotation
            if ( hist.count( hname ) > 0 ){
                ((TH2*)hist[ hname ]) -> Fill( sx, sy, v );

                // HORIZONTAL Strips
                hname = TString::Format( "hDIG1L%dhG%d", q, xG ).Data();
                ((TH2*)hist[ hname ]) -> Fill( sx, sy, v );
                    // local position
                hname = TString::Format( "hDIG1pL%dhG%d", q, xG ).Data();
                ((TH2*)hist[ hname ]) -> Fill( bcx, bcy, v );
                    // generator level
                hname = TString::Format( "hGEN1pL%dhG%d", q, xG ).Data();
                ((TH2*)hist[ hname ]) -> Fill( x, y, v );

                // VERTICAL Strips
                hname = TString::Format( "hDIG1L%dvG%d", q, yG ).Data();
                ((TH2*)hist[ hname ]) -> Fill( sx, sy, v );
                    // local position
                hname = TString::Format( "hDIG1pL%dvG%d", q, yG ).Data();
                ((TH2*)hist[ hname ]) -> Fill( bcx, bcy, v );
                    // generator level
                hname = TString::Format( "hGEN1pL%dvG%d", q, yG ).Data();
                ((TH2*)hist[ hname ]) -> Fill( x, y, v );
            }

        }
    }
}

void saturate( TH1 * h ){
    for ( int i = 1; i < h->GetNbinsX(); i++ ){
        if ( h->GetBinContent( i ) < sat_above ) continue;

        h->SetBinContent( i, (int)(sat_above - 1) );
    }
}

void make_1d( string base, int iMod, int iG, string md="" ){

    // HORIZONTAL
    string n = TString::Format( "%sL%dvG%d", md.c_str(), iMod, iG ).Data();

    hist[ (base + "2" + n) ] -> Add( hist[(base + "0" + n)] );
    hist[ (base + "2" + n) ] -> Add( hist[(base + "1" + n)] );
            // Project to 1D
    hist[ (base + "3" + n) ] = ((TH2*)hist[ (base + "2" + n) ])->ProjectionY( (base + "3" + n).c_str() );

    // VERTICAL
    n = TString::Format( "%sL%dvG%d", md.c_str(), iMod, iG ).Data();

    hist[ (base + "2" + n) ] -> Add( hist[(base + "0" + n)] );
    hist[ (base + "2" + n) ] -> Add( hist[(base + "1" + n)] );
            // Project to 1D
    hist[ (base + "3" + n) ] = ((TH2*)hist[ (base + "2" + n) ])->ProjectionX( (base + "3" + n).c_str() );
            
            
            // apply saturation on final ADC (maybe not a good approximation)
    // hist[ (base + "4" + n) ] = (TH1*)hist[ (base + "3" + n) ]->Clone( (base + "4" + n).c_str() );
    // saturate( hist[ ("h" + base + "4" + n) ] );
}


int main( int argc, char** argv ){

    loguru::init(argc, argv);

    string output_name = "output.root";
    string log_name = "everything.log";

    if (argc >= 2) 
        output_name = string( argv[1] );
    if (argc >= 3) 
        log_name = string( argv[2] );
    
    // loguru::g_stderr_verbosity = 8;
    loguru::add_file( log_name.c_str() , loguru::Truncate, loguru::Verbosity_MAX);

    LOG_F( INFO, "Writing: ROOT: %s, LOG: %s", output_name.c_str(), log_name.c_str() );

    LOG_SCOPE_FUNCTION( INFO );

    gRandom = new TRandom3();

    
    hist["hquad1"] = new TH2F( "hquad1", ";x;y", 600, 0, 600, 600, 0, 600 );
    hist["hstgc"] = new TH2F( "hstgc", ";x;y", 1400, -700, 700, 1200, -600, 600 );

    hist["hGEN0"] = new TH2F( "hGEN0", ";x;y", 1400, -700, 700, 1200, -600, 600 );
    hist["hGEN1"] = new TH2F( "hGEN1", ";x;y", 1400, -700, 700, 1200, -600, 600 );
    hist["hGEN2"] = new TH2F( "hGEN2", ";x;y", 1400, -700, 700, 1200, -600, 600 );
    hist["hDIG0"] = new TH2F( "hDIG0", ";strip x; strip y", 400, -200, 200, 400, -200, 200 );
    hist["hDIG1"] = new TH2F( "hDIG1", ";strip x; strip y", 400, -200, 200, 400, -200, 200 );
    hist["hDIG2"] = new TH2F( "hDIG2", ";strip x; strip y", 400, -200, 200, 400, -200, 200 );
    for ( int i = 0; i < 3; i++ ){
        for ( int iMod = 1; iMod < 5; iMod++ ){
            string n = TString::Format( "hDIG%dL%d", i, iMod ).Data();
            LOG_F( INFO, "Making : %s", n.c_str() );
            hist[n] = new TH2F( n.c_str(), TString::Format("Module %d;local strip x; local strip y", iMod), 200, 0, 200, 200, 0, 200 );
            for ( int ixG = 0; ixG < 3; ixG++ ){
                string n = TString::Format( "hDIG%dL%dhG%d", i, iMod, ixG ).Data();
                hist[n] = new TH2F( n.c_str(), TString::Format("Module %d strip_group %d;local strip x; local strip y", iMod, ixG), 200, 0, 200, 200, 0, 200 );

                n = TString::Format( "hDIG%dpL%dhG%d", i, iMod, ixG ).Data();
                hist[n] = new TH2F( n.c_str(), TString::Format("Module %d strip_group %d;local x; local y", iMod, ixG), 170, 0, strip_pitch * 170, 170, 0, strip_pitch * 170 );

                // GENERATOR LEVEL
                n = TString::Format( "hGEN%dpL%dhG%d", i, iMod, ixG ).Data();
                hist[n] = new TH2F( n.c_str(), TString::Format("Module %d strip_group %d;local strip x; local strip y", iMod, ixG), 170, 0, strip_pitch * 170, 170, 0, strip_pitch * 170 );
                
            }
            for ( int iyG = 0; iyG < 3; iyG++ ){
                string n = TString::Format( "hDIG%dL%dvG%d", i, iMod, iyG ).Data();
                hist[n] = new TH2F( n.c_str(), TString::Format("Module %d yGroup %d;local strip x; local strip y", iMod, iyG), 200, 0, 200, 200, 0, 200 );

                n = TString::Format( "hDIG%dpL%dvG%d", i, iMod, iyG ).Data();
                hist[n] = new TH2F( n.c_str(), TString::Format("Module %d yGroup %d;local x; local y", iMod, iyG), 170, 0, strip_pitch * 170, 170, 0, strip_pitch * 170 );

                // GENERATOR LEVEL
                n = TString::Format( "hGEN%dpL%dvG%d", i, iMod, iyG ).Data();
                hist[n] = new TH2F( n.c_str(), TString::Format("Module %d strip_group %d;local strip x; local strip y", iMod, iyG), 170, 0, strip_pitch * 170, 170, 0, strip_pitch * 170 );
                
            }
        }
    }

    fClusterProfile = new TF1( "fClusterProfile", "gaus" );
    fClusterProfile->SetParameters( 1.0, 0, 1.4 * 3.2 );
    fClusterProfile->SetRange( -500, 500 );


    // Just show that the bounds mapping is working
    for ( int i = 0; i <1000000; i++ ){

        float rx = gRandom->Rndm() * 1200 - 600;
        float ry = gRandom->Rndm() * 1200 - 600;

        if ( in_bounds( rx, ry ) )
            hist["hquad1"]->Fill( rx, ry );
        if ( in_bounds_quad( rx, ry ) > 0 ){
            hist["hstgc"]->Fill( rx, ry );

            // ****************************************************************
            // fill noise
            // ****************************************************************
            if ( gRandom->Rndm() < noise_prob ){
                float nv = gRandom->Rndm() * noise_level;
                ((TH2*)hist["hGEN2"])->Fill( rx, ry,  nv );
                int sx = -1, sy = -1;
                map_to_strip( rx, ry, sx, sy );
                ((TH2*)hist["hDIG0"])->Fill( sx, sy, nv );

                int q = map_to_local_strip( rx, ry, sx, sy );
                string hname = TString::Format( "hDIG0L%d", q ).Data();

                int xG = strip_group( sx, sy );
                int yG = strip_group( sy, sx ); // not a typo, just reuse function with rotation

                if ( hist.count( hname ) > 0 ){
                    ((TH2*)hist[ hname ]) -> Fill( sx, sy, nv );

                    hname = TString::Format( "hDIG0L%dhG%d", q, xG ).Data();
                    ((TH2*)hist[ hname ]) -> Fill( sx, sy, nv );

                    hname = TString::Format( "hDIG0pL%dhG%d", q, xG ).Data();
                    ((TH2*)hist[ hname ]) -> Fill( rx, ry, nv );

                    hname = TString::Format( "hDIG0L%dvG%d", q, yG ).Data();
                    ((TH2*)hist[ hname ]) -> Fill( sx, sy, nv );

                    hname = TString::Format( "hDIG0pL%dvG%d", q, yG ).Data();
                    ((TH2*)hist[ hname ]) -> Fill( rx, ry, nv );

                }
            }
        }

    }


    int nClusters = 0;
    while( nClusters < n_clusters_to_gen ){
        float rx = gRandom->Rndm() * 1200 - 600;
        float ry = gRandom->Rndm() * 1200 - 600;
        if ( in_bounds_quad( rx, ry ) > 0 ){
            ((TH2*)hist["hGEN0"])->Fill( rx, ry, 1 );
            fill_cluster_GEN( rx, ry );
            hit_graph->SetPoint(nClusters, rx, ry);
            nClusters++;
        }
    }

    hist["hGEN2"]->Add( hist["hGEN1"] );

    //  hDIG2 = hDIG0 + hDIG1 (noise + signal)
    hist["hDIG2"]->Add( hist[ "hDIG0" ] );
    hist["hDIG2"]->Add( hist[ "hDIG1" ] );

    // add for each module and strip group
    // also project them to make the 1D
    for ( int iMod = 1; iMod <= 4; iMod ++ ){
        string n = TString::Format( "L%d", iMod ).Data();
        hist[ ("hDIG2" + n) ] -> Add( hist[("hDIG0" + n)] );
        hist[ ("hDIG2" + n) ] -> Add( hist[("hDIG1" + n)] );

        for ( int iG = 0; iG < 3; iG++ ){
            // HORIZONTAL 
                // local position
            // string n = TString::Format( "pL%dhG%d", iMod, iG ).Data();
            make_1d( "hDIG", iMod, iG );
            make_1d( "hDIG", iMod, iG, "p" );

            make_1d( "hGEN", iMod, iG, "p" );
        }
    }


    TFile * fout = new TFile( output_name.c_str(), "RECREATE" );
    fout->cd();

    for ( auto nh : hist ){
        if ( nullptr != nh.second )
            nh.second->Write();
    }


    fClusterProfile->Write();
    hit_graph->SetName( "hitGraph" );
    hit_graph->Write();
    fout->Write();
    fout->Close();

}