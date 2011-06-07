#!/usr/bin/env python

#---------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV
##


#---------------------------------------------------------------------------
class solver( object ):
    def __init__( self, runTime, U, p, phi, transportProperties ):
        self.runTime_ = runTime
        self.U_ = U
        self.p_ = p
        self.phi_ = phi
        self.transportProperties_ = transportProperties
        self.mesh_ = U.mesh()
        self.pRefCell_ = 0
        self.pRefValue_ = 0.0
        self.pressureRes_ = 0.0
        self.velocityRes_ = 0.0
        
        from Foam.finiteVolume import setRefCell
        from Foam.OpenFOAM import word
        self.pRefCell_, self.pRefValue_ = setRefCell( self.p_, 
                                                      self.mesh_.solutionDict().subDict( word( "PISO" ) ),
                                                      self.pRefCell_,
                                                      self.pRefValue_ )

        pass
    
    
    #-----------------------------------------------------------------------
    def step( self ):
        from Foam.OpenFOAM import ext_Info, nl
        if self.runTime_.end():
            ext_Info() << "Reached end time.  Exiting" << nl
            return
        
        self.runTime_.increment()
                
        # Read transport properties
        from Foam.OpenFOAM import dimensionedScalar, word
        nu = dimensionedScalar( self.transportProperties_.lookup( word( "nu" ) ) )
                
        if self.mesh_.nInternalFaces():
            SfUfbyDelta = self.mesh_.deltaCoeffs()*self.phi_.mag()
            CoNum = ( SfUfbyDelta / self.mesh_.magSf() ).ext_max().value() * self.runTime_.deltaT().value()
            meanCoNum = ( SfUfbyDelta.sum() / self.mesh_.magSf().sum() ).value() * self.runTime_.deltaT().value()
            
            ext_Info() << "Courant Number mean: " << meanCoNum << " max: " << CoNum << nl
            pass

        # Read controls
        piso = self.mesh_.solutionDict().subDict( word( "PISO" ) )
        from Foam.OpenFOAM import readInt        
        nCorr = readInt( piso.lookup( word( "nCorrectors" ) ) )
                
        nNonOrthCorr = 0
        if piso.found( word( "nNonOrthogonalCorrectors" ) ):
            nNonOrthCorr = readInt( piso.lookup( word( "nNonOrthogonalCorrectors" ) ) )
            pass
           
        from Foam import fvm, fvc
        UEqn = fvm.ddt( self.U_ ) + fvm.div( self.phi_, self.U_ ) - fvm.laplacian( nu, self.U_ )
        
        from Foam.finiteVolume import solve
        self.velocityRes_ = solve( UEqn == -fvc.grad( self.p_ ) ).initialResidual()
                
        # --- PISO loop
                
        for corr in range( nCorr ):
            rUA = 1.0 / UEqn.A()
                    
            self.U_.ext_assign( rUA * UEqn.H() )
            self.phi_.ext_assign( ( fvc.interpolate( self.U_ ) & self.mesh_.Sf() ) + fvc.ddtPhiCorr( rUA, self.U_, self.phi_) )
                    
#           adjustPhi(phi_, U_, p_);
                    
            for nonOrth in range( nNonOrthCorr +1 ):
                pEqn =  fvm.laplacian( rUA, self.p_ ) == fvc.div( self.phi_ )
                
                pEqn.setReference( self.pRefCell_, self.pRefValue_);
                self.pressureRes_ = pEqn.solve().initialResidual()
                
                if nonOrth == nNonOrthCorr:
                   self.phi_.ext_assign( self.phi_ - pEqn.flux() )
                   pass
                pass 
                    
            # Continuity errors
            contErr = fvc.div( self.phi_ )
                    
            sumLocalContErr = self.runTime_.deltaT().value() * contErr.mag().weightedAverage( self.mesh_.V() ).value() 
                   
            globalContErr = self.runTime_.deltaT().value() * contErr.weightedAverage( self.mesh_.V() ).value()
                    
            ext_Info() << "time step continuity errors : sum local = " << sumLocalContErr << ", global = " << globalContErr << nl
                    
            # Correct velocity
            self.U_.ext_assign( self.U_ - rUA * fvc.grad( self.p_ ) )
            self.U_.correctBoundaryConditions()
             
            pass


    #-----------------------------------------------------------------------
    #- Residuals
    def pressureRes( self ):
        return self.pressureRes_
            
    def velocityRes( self ):
        return self.velocityRes_;


#---------------------------------------------------------------------------
