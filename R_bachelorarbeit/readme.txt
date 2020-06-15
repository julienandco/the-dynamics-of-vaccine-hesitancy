##readme für die Datei eigenvalue_diagram_generation.R:

Wahl von my_fav_param:

#alpha gibt nur reelle EV
#bbeta gibt nur reelle EV
#c bleibt sehr nah an reeller achse (sieht aus wie alpha und beta, hat nur einen im.wert bei -1.6+- 0.6i ugf)
#a gibt nur reelle EV
#für B bekomme ich schönes blatt, aber trotzdem kein iomega...
#n1 hat ab Re >= -1.0 nur noch reelle EV
#n2 im bereich der y-achse nur noch reelle EV
#theta1 konvergiert nicht mehr ab tt = 1.74, davor nur reell
#theta2 konvergiert nicht mehr ab tt = 1.22, davor nur reell


Wählen also B, da schönes Blatt. Wenn wir aber andere Parameter zu Beginn ändern und danach mit B iterieren:

#default parameter:

bbeta = 5
alpha = 0.45
B = 2
a = 0.1
c = 0.3

theta1 = 1
theta2 = theta1
n1 = 10
n2 = 10


###mit zufällig gewählten Parametern:

#großes a gibt doppelblatt (kleineres Blatt im größeren integriert)
#großes a und großes c gibt einen Pfeil
#grßes a, c und nx gibt nur reelle EW
#großes n1 gibt nur reelle EW
#kleines n1 streckt Blatt
#großes n2 streckt Blatt
#kleines n2 gibt nur relle EW
#großes bbeta streckt blatt
#kleines bbeta staucht blatt
# 0 <= alpha <= 0.1 gibt winziges blatt auf reeller Achse
# 0.15 <= alpha <= 0.45 gibt ganz verschiedene Objekte, die sich mehr und mehr dem Blatt nähern: 0.15 = zwei Ovale auf reeller Achse, 0.2 = Raumschiff, ab 0.25 bis 0.45 nähern wir uns dem Blatt
# 0.45 < alpha <= 1 gibt aufgeblähtes Blatt
# 1 < alpha: Blatt schrumpft langsam, bis es nur noch reelle EW gibt




###mit kalibrierten Parametern (sodass x=1/2 stat. Pkt.):

#großes a konvergiert nicht...
#großes a und großes c konvergiert nicht...
#alpha zeigt ähnliches Verhalten wie oben: 0.1= zwei Blätter auf reeller Achse , 0.15 = Birne, 0.2 = Birne uvm.


