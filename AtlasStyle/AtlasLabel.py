from ROOT import TLatex, gPad

def ATLASLabel(x,y,text,color=1):

      l = TLatex()
      l.SetNDC()
      l.SetTextFont(72)
      l.SetTextColor(color)
      #l.SetTextSize(0.06)
      l.SetTextSize(0.04)

      delx = 0.115*550*gPad.GetWh()/(472*gPad.GetWw())
      
      #delx = 0.115*500*gPad.GetWh()/(472*gPad.GetWw())

      l.DrawLatex(x,y,"ATLAS")

      if (text):
          p = TLatex()
          p.SetNDC()
          #p.SetTextSize(0.06)
          p.SetTextSize(0.04)
          p.SetTextFont(42)
          p.SetTextColor(color)
          p.DrawLatex(x+delx,y,text)

def myText(x,y,text, color = 1):
    l = TLatex()
    l.SetTextSize(0.04)
    #l.SetTextSize(0.025)
    l.SetNDC()
    l.SetTextColor(color)
    l.DrawLatex(x,y,text)
