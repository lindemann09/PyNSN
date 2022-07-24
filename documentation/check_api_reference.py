import sys, os
import inspect
import json

sys.path.insert(0, os.path.abspath('..'))
import pynsn
from pynsn.image import pil_image,mpl_figure

print("PyNSN Version: " + pynsn.__version__)


def inspect_all_members(item, rtn_dict):

    def inspect_members(sub_item, rtn_dict):
        members = inspect.getmembers(eval(sub_item))
        modules = []
        for member in members:
            if member[0][0:1] != '_':
                if inspect.ismodule(member[1]):
                    modules.append(sub_item + "." + member[0])
                elif inspect.isclass(member[1]):
                    rtn_dict["classes"].append(sub_item + "." + member[0])
                elif inspect.isfunction(member[1]):
                    rtn_dict["functions"].append(sub_item + "." + member[0])
                else:
                    rtn_dict["attributes"].append(sub_item + "." + member[0])

        for m in modules:
            rtn_dict["modules"].append(m)
            rtn_dict = inspect_members(m, rtn_dict) # recursive

        return rtn_dict

    return inspect_members(item, rtn_dict)


def get_undocumented_members():
    # get all members
    members = {"modules": [], "classes": [], "functions": [], "attributes": []}
    members = inspect_all_members("pynsn", members)
    del members["modules"]
    # del members["attributes"] # todo

    # get autodoc: all lines with ..autodoc
    autodoc_txt = ""
    for flname in os.listdir("."):
        if flname.endswith(".rst"):
            with open(flname, "r") as fl:
                for l in fl.readlines():
                    if l.find(".. auto") >= 0:
                        autodoc_txt += l.strip() + " "

    undoc = []
    for k in members.keys():
        for m in members[k]:
            if autodoc_txt.find(m) < 0:
                print("!!Warning!! Undocumented library member: "+ m)
                undoc.append(m)
    return undoc

get_undocumented_members()
