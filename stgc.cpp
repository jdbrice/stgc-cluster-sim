

#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom3.h"

#include <iostream>
#include <map>

#define LOGURU_IMPLEMENTATION 1
#include "loguru.h"

// in mm
const float pent_base = 537.0;
const float pent_nib = 179.0;
const float pent_shift = 101.6;
const int n_clusters_to_gen = 200;
const float noise_prob = 0.7; // 0 -1, 1 is max noise
const float noise_level = 10; // ADC
const float strip_pitch = 3.2;
const size_t sat_above = 1024; // ADC
const float cluster_max_adc = 300; // ADC

TF1 * fClusterProfile = nullptr;

using namespace std;
map<string, TH1*> hist;

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
    } else if ( in_bounds( x - pent_shift, -y ) ){
        return 4;
    }

    return 0;
}

void map_to_strip( float x, float y, int &sx, int &sy ){
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
    
}

void fill_cluster_GEN( float x, float y, TH2 * hGEN, TH2 * hDIG ){
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
            // if ( q == 1 )
            hDIG->Fill( sx, sy, (int)v );

        }
    }
}

int main( int argc, char** argv ){

    loguru::init(argc, argv);
    // loguru::g_stderr_verbosity = 8;
    loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);

    LOG_SCOPE_FUNCTION( INFO );

    gRandom = new TRandom3();

    

    hist["hquad1"] = new TH2F( "hquad1", ";x;y", 600, 0, 600, 600, 0, 600 );
    hist["hstgc"] = new TH2F( "hstgc", ";x;y", 1400, -700, 700, 1200, -600, 600 );

    hist["hGEN0"] = new TH2F( "hGEN0", ";x;y", 1400, -700, 700, 1200, -600, 600 );
    hist["hGEN1"] = new TH2F( "hGEN1", ";x;y", 1400, -700, 700, 1200, -600, 600 );
    hist["hGEN2"] = new TH2F( "hGEN2", ";x;y", 1400, -700, 700, 1200, -600, 600 );
    hist["hDIG"] = new TH2F( "hDIG", ";x;y", 400, -200, 200, 400, -200, 200 );

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


            // fill noise
            if ( gRandom->Rndm() < noise_prob ){
                float nv = gRandom->Rndm() * noise_level;
                ((TH2*)hist["hGEN2"])->Fill( rx, ry,  nv );
                int sx = -1, sy = -1;
                map_to_strip( rx, ry, sx, sy );
                ((TH2*)hist["hDIG"])->Fill( sx, sy, nv );
            }

            
        }

    }


    int nClusters = 0;
    
    while( nClusters < n_clusters_to_gen ){
        float rx = gRandom->Rndm() * 1200 - 600;
        float ry = gRandom->Rndm() * 1200 - 600;
        if ( in_bounds_quad( rx, ry ) > 0 ){
            ((TH2*)hist["hGEN0"])->Fill( rx, ry, 1 );
            fill_cluster_GEN( rx, ry, (TH2*)hist["hGEN1"], (TH2*)hist["hDIG"] );
            nClusters++;
        }
    }

    hist["hGEN2"]->Add( hist["hGEN1"] );

    // check all the bins of hhDIG and apply saturation
    for ( int ix = 1; ix < hist["hDIG"]->GetXaxis()->GetNbins(); ix++ ){
        for ( int iy = 1; iy < ((TH2*)hist["hDIG"])->GetYaxis()->GetNbins(); iy++ ){
            float bv = ((TH2*)hist["hDIG"])->GetBinContent( ix, iy );
            if ( bv >= sat_above ){
                ((TH2*)hist["hDIG"])->SetBinContent( ix, iy, sat_above );
                ((TH2*)hist["hDIG"])->SetBinError( ix, iy, 0.001 );
            }
        }
    }


    TFile * fout = new TFile( "output.root", "RECREATE" );
    fout->cd();

    for ( auto nh : hist ){
        if ( nullptr != nh.second )
            nh.second->Write();
    }


    fClusterProfile->Write();

    fout->Write();
    fout->Close();

}