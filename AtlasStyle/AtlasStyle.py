# $Id: AtlasStyle.py 200750 2011-08-19 09:50:04Z skreiss $
#
# Module implementing a TStyle that follows the ATLAS recommendations.
#

# PyROOT is needed:
import ROOT

##
# @short An ATLAS-specific style
#
# This is yet another implementation for the common ATLAS style. The user
# basically just needs to import this module, and from there on the plots
# produced will follow the atlas plot style.
#
class AtlasStyle( ROOT.TStyle ):

    ##
    # @short Object constructor
    #
    # The constructor just initializes the underlying TStyle object, and
    # calls the configure() function to set up the ATLAS style.
    #
    # The parameters of the constructor should just be ignored in 99.9%
    # of the cases.
    #
    # @param name The name given to the style
    # @param title The title given to the style
    def __init__( self, name = "AtlasStyle", title = "ATLAS style object" ):

        print "Setting ATLAS style"

        # Initialise the base class:
        ROOT.TStyle.__init__( self, name, title )
        self.SetName( name )
        self.SetTitle( title )

        # Call the configure function for setting up the style:
        self.configure()

        return

    ##
    # @short Configure the object for the ATLAS style
    #
    # This function actually takes care of setting up the ATLAS style.
    # I implemented it based on a C++ TStyle object, which in turn was
    # implemented based on a central piece of CINT macro.
    def configure( self ):

        # Tell the user what we're doing:
        self.Info( "configure", "Configuring default ATLAS style" )

        # Use plain black on white colors:
        icol = 0
        self.SetFrameBorderMode( 0 )
        self.SetFrameFillColor( icol )
        self.SetFrameFillStyle( 0 )
        self.SetCanvasBorderMode( 0 )
        self.SetPadBorderMode( 0 )
        self.SetPadColor( icol )
        self.SetCanvasColor( icol )
        self.SetStatColor( icol )

        # Set the paper and margin sizes:
        self.SetPaperSize( 20, 26 )
        self.SetPadTopMargin( 0.07 )
        self.SetPadRightMargin( 0.05 )
        self.SetPadBottomMargin( 0.16 )
        self.SetPadLeftMargin( 0.16 )

        # set title offsets (for axis label)
        self.SetTitleXOffset(1.3);
        self.SetTitleYOffset(1.3);

        # Use large fonts:
        font_type = 42
        font_size = 0.05
        self.SetTextFont( font_type )
        self.SetTextSize( font_size )
        self.SetLegendFont( font_type )
        #self.SetLegendSize( font_size )
        self.SetLabelFont( font_type, "x" )
        self.SetLabelSize( font_size, "x" )
        self.SetTitleFont( font_type, "x" )
        self.SetTitleSize( font_size, "x" )
        self.SetLabelFont( font_type, "y" )
        self.SetLabelSize( font_size, "y" )
        self.SetTitleFont( font_type, "y" )
        self.SetTitleSize( font_size, "y" )
        self.SetLabelFont( font_type, "z" )
        self.SetLabelSize( font_size, "z" )
        self.SetTitleFont( font_type, "z" )
        self.SetTitleSize( font_size, "z" )

        # Use bold lines and markers:
        self.SetMarkerStyle( 20 )
        self.SetMarkerSize( 1.2 )
        self.SetHistLineWidth( 2 )
        self.SetLineStyleString( 2, "[12 12]" )

        # Do not display any of the standard histogram decorations:
        self.SetOptTitle( 0 )
        self.SetOptStat( 0 )
        self.SetOptFit( 0 )

        # Put tick marks on top and rhs of the plots:
        self.SetPadTickX( 1 )
        self.SetPadTickY( 1 )

        return

# Tell ROOT to use this style:
style = AtlasStyle()
ROOT.gROOT.SetStyle( style.GetName() )
# Force style applies styles to histograms made under a different style regim
ROOT.gROOT.ForceStyle()
ROOT.TGaxis.SetMaxDigits( 4 )
#ROOT.gStyle.SetPalette(1)
