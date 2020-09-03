from typing import List

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

def plot(forward: List[float], reverse: List[float], font: str, output: str):
    f = matplotlib.font_manager.FontProperties(fname = font)
    matplotlib.pyplot.xlabel("distance from footprint (bp)", fontproperties = f)
    matplotlib.pyplot.ylabel("average bias-corrected Tn5 insertions", fontproperties = f)
    matplotlib.pyplot.gca().spines['right'].set_visible(False)
    matplotlib.pyplot.gca().spines['top'].set_visible(False)
    matplotlib.pyplot.plot(range(-500, 500), forward, label = "foward strand")
    matplotlib.pyplot.plot(range(-500, 500), reverse, label = "reverse strand")
    matplotlib.pyplot.xticks([ -500, -250, 0, 250, 500 ], fontproperties = f)
    matplotlib.pyplot.yticks(range(0, 1000, 100), fontproperties = f)
    matplotlib.pyplot.legend(prop = f)
    matplotlib.pyplot.savefig(output)
