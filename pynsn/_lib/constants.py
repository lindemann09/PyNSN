from __future__ import unicode_literals

P_AREA = u"Total surface area"
P_DIAMETER = u"Mean dot diameter"
P_CIRCUMFERENCE = u"Total circumference"
P_CONVEX_HULL = u"Convex hull"
P_DENSITY = u"Density"
PROPERTY_DEPENDENCIES = [[P_CONVEX_HULL, P_DENSITY], [P_DIAMETER, P_CIRCUMFERENCE, P_AREA]]


