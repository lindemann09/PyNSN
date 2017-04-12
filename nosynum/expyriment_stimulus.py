from __future__ import absolute_import, print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from expyriment.stimuli import Canvas, Circle, Line, Picture

def create_stimulus(dot_array, area_colour=None,
                    convex_hull_colour=None,
                    antialiasing=None,  #TODO
                               background_stimulus_expyriment=None):
    if background_stimulus_expyriment is None:
        canvas = Canvas(size=(dot_array._stimulus_area_radius * 2,) * 2)
    else:
        canvas = background_stimulus_expyriment.copy()

    if area_colour is not None:
        Circle(radius=dot_array._stimulus_area_radius,
               colour=area_colour).plot(canvas)
    if convex_hull_colour is not None:
        # plot convey hull
        hull = dot_array.convex_hull
        hull.append(hull[0])
        last = None
        for p in hull:
            if last is not None:
                Line(start_point=last, end_point=p, line_width=2,
                     colour=convex_hull_colour).plot(canvas)
            last = p

    # plot dots
    for d in dot_array.dots:
        if d.picture is not None:
            Picture(filename=d.picture,
                    position=d.xy).plot(canvas)
        else:
            Circle(radius=d.radius, colour=d.colour,
                   line_width=0, position=d.xy).plot(canvas)

    return canvas