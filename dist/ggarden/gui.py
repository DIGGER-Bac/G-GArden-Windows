import datetime
import itertools
import csv
import os

import numpy as np
import pandas as pd
from kivy.app import App
from kivy.uix.textinput import TextInput
from kivy.properties import DictProperty, ObjectProperty, ListProperty, StringProperty, NumericProperty
from kivy.uix.floatlayout import FloatLayout
from kivy.uix.togglebutton import ToggleButton
from kivy.uix.widget import Widget
from kivy.uix.button import Button
from kivy.uix.dropdown import DropDown
from kivy.uix.popup import Popup
from kivy.uix.checkbox import CheckBox
from kivy.uix.slider import Slider
from datasets import registered_fragments, translations
from tool_main import all_in_one_better
from helpers import build_prim_coll, detect_illegal_res_sites, find_overhang, input_handler, translate, ligation_fidelity
from kivy.uix.label import Label
from read_input import read_in_file_fasta, read_seq_from_csv, check_meta_inputs
from kivy.uix.spinner import Spinner, SpinnerOption
from datasets import res_2s
import shutil
from collections import Counter
import webbrowser
import pyperclip3

# contains all selected fragments from MultiSelectSpinner instances as nested lists
selection = []

# contains all custom overhangs from OhInput instances
custom_ohs = []

# contains all custom primer processings from PrimerSpinner instances
custom_primer = []


class OptimButton(Button):

    ohs = StringProperty("")

    def optim_link(self):
        # pyclip.copy("ACCA,GTTG,AACG,TTGA,AGTC")
        webbrowser.open(url="https://ggtools.neb.com/getset/run.cgi", new=0, autoraise=True)


class OptimizeCheck(CheckBox):
    def __init__(self, **kwargs):
        super(OptimizeCheck, self).__init__(**kwargs)
        self.bind(active=self.on_checkbox_active)

    def on_checkbox_active(self, instance, value):
        if value:
            print('The checkbox', instance.text, 'is active')
        else:
            print('The checkbox', instance.text, 'is inactive')


class MyDropDown(DropDown):
    pass


class ProtocolButton(Button):
    file_conns = DictProperty({"T4 01h 25C": "FileS01_T4_01h_25C.xlsx",
                               "T4 01h 37C": "FileS02_T4_01h_37C.xlsx",
                               "T4 18h 25C": "FileS03_T4_18h_25C.xlsx",
                               "T4 18h 37C": "FileS04_T4_18h_37C.xlsx"})

    def create_products(self, selected):
        products = list(itertools.product(*selected))
        uniques = itertools.chain.from_iterable(selected)
        for fragment in uniques:
            check_seq = registered_fragments[fragment][2]
            check_seq = check_seq.replace("[", "")
            check_seq = check_seq.replace("]", "")
            if detect_illegal_res_sites(querry=check_seq, motif=res_2s[self.parent.ids["meta"].meta_dict["Restriction Enzyme"]][0]):
                raise ValueError(f"Fragment {fragment} contains recognition site for type IIS enzyme")
        return products

    def visualize(self):
        try:
            app = App.get_running_app()
            length = int(app.root.ids["meta"].ids["oh_len"].ids["txt"].text)
            res_enz = app.root.ids["meta"].ids["res_enz"].ids["txt"].text
            print(custom_ohs)
            oh_collection = []
            products = self.create_products(selected=selection)
            print(products)
            for product in products:
                prod_ohs = []
                for idx, fragment in enumerate(product):
                    if len(custom_ohs[idx]) == 0:
                        if idx + 1 < len(product):
                            conn = f"{registered_fragments[fragment][1][:4]}-{registered_fragments[product[idx+1]][1][:4]}"
                            seqs = [["placeholder", registered_fragments[fragment][2]], ["placeholder", registered_fragments[product[idx+1]][2]]]
                            rel_seq = input_handler(conn=conn, idx=0, content=seqs)
                            prod_ohs.append(find_overhang(seq=rel_seq, conn_type=conn, length=length, res_enz=res_enz))
                        else:
                            conn = f"{registered_fragments[fragment][1][:4]}-{registered_fragments[product[0]][1][:4]}"
                            seqs = [["placeholder", registered_fragments[fragment][2]], ["placeholder", registered_fragments[product[0]][2]]]
                            rel_seq = input_handler(conn=conn, idx=0, content=seqs)
                            prod_ohs.append(find_overhang(seq=rel_seq, conn_type=conn, length=length, res_enz=res_enz))
                    else:
                        prod_ohs.append([translate(custom_ohs[idx]), custom_ohs[idx]])
                oh_collection.append(prod_ohs)
            # excel_data = pd.read_excel(f'sb8b00333_si_002/Supplemental Data/{self.file_conns[self.text]}', header=0, index_col=0)
            # data = pd.DataFrame(excel_data)
            # for oh_set in oh_collection:
            #     unlisted_oh_set = list(itertools.chain(*oh_set))
            #     print(unlisted_oh_set)
            #     fids, fid = ligation_fidelity(ohs=unlisted_oh_set, mat=data)
            #     print(fid)
            oh_view = OhViewPopup(ohs=oh_collection, prods=products, mat_file="FileS03_T4_18h_25C.xlsx")
            oh_view.open()
        except Exception as exc:
            self.parent.ids["vis_out"].text = f"[ref=reset]{str(exc)}[/ref]"


class ProtocolSpinner(Spinner):

    file_conns = DictProperty({"T4 01h 25C": "FileS01_T4_01h_25C.xlsx",
                               "T4 01h 37C": "FileS02_T4_01h_37C.xlsx",
                               "T4 18h 25C": "FileS03_T4_18h_25C.xlsx",
                               "T4 18h 37C": "FileS04_T4_18h_37C.xlsx"})

    def create_products(self, selected):
        products = list(itertools.product(*selected))
        uniques = itertools.chain.from_iterable(selected)
        for fragment in uniques:
            check_seq = registered_fragments[fragment][2]
            check_seq = check_seq.replace("[", "")
            check_seq = check_seq.replace("]", "")
            if detect_illegal_res_sites(querry=check_seq, motif=res_2s[self.parent.ids["meta"].meta_dict["Restriction Enzyme"]][0]):
                raise ValueError(f"Fragment {fragment} contains recognition site for type IIS enzyme")
        return products

    def visualize(self):
        try:
            if self.text == "Scarless" or self.text in self.file_conns:
                print(custom_ohs)
                oh_collection = []
                products = self.create_products(selected=selection)
                print(products)
                for product in products:
                    prod_ohs = []
                    for idx, fragment in enumerate(product):
                        if len(custom_ohs[idx]) == 0:
                            if idx + 1 < len(product):
                                conn = f"{registered_fragments[fragment][1][:4]}-{registered_fragments[product[idx+1]][1][:4]}"
                                seqs = [["placeholder", registered_fragments[fragment][2]], ["placeholder", registered_fragments[product[idx+1]][2]]]
                                rel_seq = input_handler(conn=conn, idx=0, content=seqs)
                                prod_ohs.append(find_overhang(seq=rel_seq, conn_type=conn, length=4, res_enz="BbsI"))
                            else:
                                conn = f"{registered_fragments[fragment][1][:4]}-{registered_fragments[product[0]][1][:4]}"
                                seqs = [["placeholder", registered_fragments[fragment][2]], ["placeholder", registered_fragments[product[0]][2]]]
                                rel_seq = input_handler(conn=conn, idx=0, content=seqs)
                                prod_ohs.append(find_overhang(seq=rel_seq, conn_type=conn, length=4, res_enz="BbsI"))
                        else:
                            prod_ohs.append([translate(custom_ohs[idx]), custom_ohs[idx]])
                    oh_collection.append(prod_ohs)
                if self.text in self.file_conns:
                    oh_view = OhViewPopup(ohs=oh_collection, prods=products, mat_file=self.file_conns[self.text])
                else:
                    oh_view = OhViewPopup(ohs=oh_collection, prods=products)
                oh_view.open()
            elif self.text == "Fidelity":
                webbrowser.open(url="https://ggtools.neb.com/getset/run.cgi", new=0, autoraise=True)
            else:
                pass
        except Exception as exc:
            self.parent.ids["vis_out"].text = f"[ref=reset]{str(exc)}[/ref]"


class FidelitySpinner(Spinner):
    pass


########################################################################################################################
# Button to delete Dropdown instances ##################################################################################
########################################################################################################################


class VisualizeOhs(Button):
    def create_products(self, selected):
        products = list(itertools.product(*selected))
        uniques = itertools.chain.from_iterable(selected)
        for fragment in uniques:
            check_seq = registered_fragments[fragment][2]
            check_seq = check_seq.replace("[", "")
            check_seq = check_seq.replace("]", "")
            if detect_illegal_res_sites(querry=check_seq, motif=res_2s[self.parent.ids["meta"].meta_dict["Restriction Enzyme"]][0]):
                raise ValueError(f"Fragment {fragment} contains recognition site for type IIS enzyme")
        return products

    def visualize(self):
        print(custom_ohs)
        oh_collection = []
        products = self.create_products(selected=selection)
        print(products)
        for product in products:
            prod_ohs = []
            for idx, fragment in enumerate(product):
                if len(custom_ohs[idx]) == 0:
                    if idx + 1 < len(product):
                        conn = f"{registered_fragments[fragment][1][:4]}-{registered_fragments[product[idx+1]][1][:4]}"
                        seqs = [["placeholder", registered_fragments[fragment][2]], ["placeholder", registered_fragments[product[idx+1]][2]]]
                        rel_seq = input_handler(conn=conn, idx=0, content=seqs)
                        prod_ohs.append(find_overhang(seq=rel_seq, conn_type=conn, length=4, res_enz="BbsI"))
                    else:
                        conn = f"{registered_fragments[fragment][1][:4]}-{registered_fragments[product[0]][1][:4]}"
                        seqs = [["placeholder", registered_fragments[fragment][2]], ["placeholder", registered_fragments[product[0]][2]]]
                        rel_seq = input_handler(conn=conn, idx=0, content=seqs)
                        prod_ohs.append(find_overhang(seq=rel_seq, conn_type=conn, length=4, res_enz="BbsI"))
                else:
                    prod_ohs.append([translate(custom_ohs[idx]), custom_ohs[idx]])
            oh_collection.append(prod_ohs)

        app = App.get_running_app()
        print(app.root.ids["meta"].ids["oh_len"].ids["txt"].text)

        excel_data = pd.read_excel(f'sb8b00333_si_002/Supplemental Data {app.root.ids["meta"].ids["oh_len"].ids["txt"].text}/FileS03_T4_18h_25C.xlsx', header=0, index_col=0)
        data = pd.DataFrame(excel_data)
        for oh_set in oh_collection:
            unlisted_oh_set = list(itertools.chain(*oh_set))
            print(unlisted_oh_set)
            fids, fid = ligation_fidelity(ohs=unlisted_oh_set, mat=data)
            print(fid)
        oh_view = OhViewPopup(ohs=oh_collection, prods=products)
        oh_view.open()


class DelDrop(Button):
    def del_last_drop(self, add_button_id, get_seq_button_id):
        max_row = [child.my_row for child in self.parent.ids["multi_layout"].children if "MultiSelectSpinner" in str(child) or "PrimerSpinner" in str(child) or "OhInput" in str(child)]
        if len(max_row) > 0:
            max_row = max(max_row)
            max_id = [child.my_id for child in self.parent.ids["multi_layout"].children if ("MultiSelectSpinner" in str(child) or "PrimerSpinner" in str(child) or "OhInput" in str(child)) and child.my_row == max_row]
            max_id = max(max_id)
            for child in self.parent.ids["multi_layout"].children[::-1]:
                if ("MultiSelectSpinner" in str(child) or "PrimerSpinner" in str(child) or "OhInput" in str(child)) and child.my_id == max_id and child.my_row == max_row:
                    self.parent.ids["multi_layout"].remove_widget(child)
                if "OhInput" in str(child) and max_id == 0 and max_row > 0 and child.my_row == max_row - 1 and child.my_id == 8:
                    self.parent.ids["multi_layout"].remove_widget(child)
                elif "OhInput" in str(child) and child.my_id == max_id - 1 and child.my_row == max_row:
                    self.parent.ids["multi_layout"].remove_widget(child)
            if self.parent.ids["multi_layout"].drop_count == 8 and self.parent.ids["multi_layout"].row_count == 2:
                self.parent.ids["multi_layout"].drop_count -= 1
                add_button_id.pos_hint = {"x": self.parent.ids["multi_layout"].dropdown_width * (self.parent.ids["multi_layout"].drop_count) + self.parent.ids["multi_layout"].xpos, "top": self.parent.ids["multi_layout"].ypos}
                get_seq_button_id.pos_hint = {"x": 0.05 + self.parent.ids["multi_layout"].dropdown_width * (self.parent.ids["multi_layout"].drop_count) + self.parent.ids["multi_layout"].xpos, "top": self.parent.ids["multi_layout"].ypos - 0.025}
                self.pos_hint = {"x": self.parent.ids["multi_layout"].dropdown_width * (self.parent.ids["multi_layout"].drop_count) + self.parent.ids["multi_layout"].xpos, "top": self.parent.ids["multi_layout"].ypos - 0.05}
                self.parent.ids["the_adder"].disabled = False
                del selection[-1]
                del custom_ohs[-1]
                del custom_primer[-1]
                print(selection, custom_ohs, custom_primer)
            elif self.parent.ids["multi_layout"].drop_count != 1 or (self.parent.ids["multi_layout"].drop_count == 1 and self.parent.ids["multi_layout"].row_count == 0):
                self.parent.ids["multi_layout"].drop_count -= 1
                add_button_id.pos_hint = {"x": self.parent.ids["multi_layout"].dropdown_width * (self.parent.ids["multi_layout"].drop_count) + self.parent.ids["multi_layout"].xpos, "top": self.parent.ids["multi_layout"].ypos}
                get_seq_button_id.pos_hint = {"x": 0.05 + self.parent.ids["multi_layout"].dropdown_width * (self.parent.ids["multi_layout"].drop_count) + self.parent.ids["multi_layout"].xpos, "top": self.parent.ids["multi_layout"].ypos - 0.025}
                self.pos_hint = {"x": self.parent.ids["multi_layout"].dropdown_width * (self.parent.ids["multi_layout"].drop_count) + self.parent.ids["multi_layout"].xpos, "top": self.parent.ids["multi_layout"].ypos - 0.05}
                del selection[-1]
                del custom_ohs[-1]
                del custom_primer[-1]
                print(selection, custom_ohs, custom_primer)
            elif self.parent.ids["multi_layout"].drop_count == 1 and self.parent.ids["multi_layout"].row_count > 0:
                self.parent.ids["multi_layout"].row_count -= 1
                self.parent.ids["multi_layout"].drop_count = 8
                self.parent.ids["multi_layout"].ypos += 3 * self.parent.ids["multi_layout"].dropdown_height
                add_button_id.pos_hint = {"x": self.parent.ids["multi_layout"].dropdown_width * (self.parent.ids["multi_layout"].drop_count) + self.parent.ids["multi_layout"].xpos, "top": self.parent.ids["multi_layout"].ypos}
                get_seq_button_id.pos_hint = {"x": 0.05 + self.parent.ids["multi_layout"].dropdown_width * (self.parent.ids["multi_layout"].drop_count) + self.parent.ids["multi_layout"].xpos, "top": self.parent.ids["multi_layout"].ypos - 0.025}
                self.pos_hint = {"x": self.parent.ids["multi_layout"].dropdown_width * (self.parent.ids["multi_layout"].drop_count) + self.parent.ids["multi_layout"].xpos, "top": self.parent.ids["multi_layout"].ypos - 0.05}
                del selection[-1]
                del custom_ohs[-1]
                del custom_primer[-1]
                print(selection, custom_ohs, custom_primer)

########################################################################################################################
# Label which displays error messages from within the main program tool_main.py ########################################
########################################################################################################################


class OutLabel(Label):
    pass

########################################################################################################################
# Button to add new dropdown instances #################################################################################
########################################################################################################################


class AddButton(Button):
    def extend_selection(self):
        selection.append([])
        custom_ohs.append("")
        custom_primer.append("")

    def disable_self(self, instance, *args):
        instance.disabled = True

########################################################################################################################
# Interactive display of outputs and output directory ##################################################################
########################################################################################################################


class LinkLabel(Label):
    pass

########################################################################################################################
# TextInput with predefined settings for usage in LabelInput ###########################################################
########################################################################################################################


class MyTextInput(TextInput):

    def __init__(self, **kwargs):
        super(MyTextInput, self).__init__(**kwargs)
        self.bind(text=self.on_text)

    def on_text(self, instance, value):
        instance.parent.parent.meta_dict[instance.parent.label] = value


########################################################################################################################
# Combination of Label and MyTextInput to get one paramemter for tool_main.py from user #####################################
########################################################################################################################


class LabelInput(Widget):

    row = NumericProperty(0)

    root_width = NumericProperty(0)

    root_height = NumericProperty(0)

    label = StringProperty()

    default = StringProperty()


class ResInput(Widget):

    row = NumericProperty(0)

    root_width = NumericProperty(0)

    root_height = NumericProperty(0)

    res = sorted(res_2s.keys())

    label = StringProperty()

    default = StringProperty()

    def change_oh_disp(self):
        self.parent.ids["oh_len"].ids["txt"].text = str(res_2s[self.ids["txt"].text][2])
        self.parent.meta_dict[self.parent.ids["oh_len"].label] = self.parent.ids["oh_len"].ids["txt"].text
        print("dsadasdasdas" + self.parent.ids["oh_len"].label)


class OhDisplay(Widget):

    row = NumericProperty(0)

    root_width = NumericProperty(0)

    root_height = NumericProperty(0)

    label = StringProperty("Overhang Length")

    default = StringProperty()

########################################################################################################################
# Combination of LabelInput instances to get all parameters for tool_main.py ################################################
########################################################################################################################


class MetaSelection(Widget):

    root_width = NumericProperty(0)

    root_height = NumericProperty(0)

    meta_dict = DictProperty({})

########################################################################################################################
# ToggleButton extension with predefined settings for usage in MultiSelectSpinner ######################################
########################################################################################################################


class MultiSelectOption(ToggleButton):
    pass

########################################################################################################################
# Spinner dropdown for primer design selection #########################################################################
########################################################################################################################


class PrimerSpinner(Spinner):

    my_id = NumericProperty()

    my_row = NumericProperty()

    def __init__(self, my_id, my_row, **kwargs):
        self.my_id = my_id
        self.my_row = my_row
        super(PrimerSpinner, self).__init__(**kwargs)
        self.bind(text=self.add_selection)

    def add_selection(self, *args):
        if self.text != "Auto":
            custom_primer[self.my_id] = self.text
        else:
            custom_primer[self.my_id] = ""

########################################################################################################################
# TextInput for custom overhang design #################################################################################
########################################################################################################################


class OhInput(TextInput):

    my_id = NumericProperty()

    my_row = NumericProperty()

    def __init__(self, my_id, my_row, **kwargs):
        self.my_id = my_id
        self.my_row = my_row
        super(OhInput, self).__init__(**kwargs)
        self.bind(text=self.add_selection)

    def add_selection(self, *args):
        custom_ohs[self.my_id] = self.text

########################################################################################################################
# Mainbutton for expansion of a fragment dropdown which is also included in this class #################################
########################################################################################################################


class MultiSelectSpinner(Button):

    dropdown = ObjectProperty(None)

    values = ListProperty([])

    selected_values = ListProperty([])

    count_selection = NumericProperty(0)

    my_id = NumericProperty()

    my_row = NumericProperty()

    def __init__(self, my_id, my_row, **kwargs):
        self.my_id = my_id
        self.my_row = my_row
        self.bind(dropdown=self.update_dropdown)
        self.bind(values=self.update_dropdown)
        self.bind(selected_values=self.add_selection)
        super(MultiSelectSpinner, self).__init__(**kwargs)
        self.bind(on_release=self.toggle_dropdown)

    def toggle_dropdown(self, *args):
        if self.dropdown.parent:
            self.dropdown.dismiss()
        else:
            self.dropdown.open(self)
            self.background_normal = 'button_up.png'

    def update_dropdown(self, *args):
        if not self.dropdown:
            self.dropdown = MyDropDown()
        values = self.values
        if values:
            if self.dropdown.children:
                self.dropdown.clear_widgets()
            for value in values:
                b = MultiSelectOption(text=value)
                b.bind(state=self.select_value)
                self.dropdown.add_widget(b)

    # managing selected values
    def select_value(self, instance, value):
        if value == 'down':
            if instance.text not in self.selected_values:
                self.selected_values.append(instance.text)
        else:
            if instance.text in self.selected_values:
                self.selected_values.remove(instance.text)

    # displaying selected values
    # instance has to stay
    def on_selected_values(self, instance, value):
        if value:
            if len(value) == 1:
                self.text = value[0]
            else:
                self.text = str(len(value))
        else:
            self.text = "0"

    def add_selection(self, instance, *args):
        print(self.my_id, self.my_row)
        own_number = self.my_id + self.my_row * 9
        selection[own_number] = self.selected_values
        frags = [registered_fragments[frag][1] for frag in self.selected_values]
        prev_id = self.my_id - 1 if self.my_id != 0 else 8
        prev_row = self.my_row if self.my_id != 0 else self.my_row - 1
        print(frags)
        # behaviour if self contains overhang restricting fragments
        if "plasmid" in frags or "dprom" in frags or "dseed" in frags or "dscaf" in frags:
            for child in self.parent.children[::-1]:
                # if oh input of previous fragment or self oh input is found: remove it
                if "OhInput" in str(child):
                    if own_number - 1 == child.my_id + child.my_row * 9 or self.my_id + self.my_row * 9 == child.my_id + child.my_row * 9:
                        self.parent.remove_widget(child)
                # if primer spinner of self is found: remove it
                elif "PrimerSpinner" in str(child) and self.my_id == child.my_id and self.my_row == child.my_row:
                    self.parent.remove_widget(child)
        # behaviour if self removes overhang restricting fragments
        else:
            # get existing primer spinner numbers and overhang numbers
            pcr_ids = [child.my_id + child.my_row * 9 for child in self.parent.children if "PrimerSpinner" in str(child)]
            oh_ids = [child.my_id + child.my_row * 9 for child in self.parent.children if "OhInput" in str(child)]
            print(oh_ids)
            if own_number not in pcr_ids:
                self.parent.add_widget(PrimerSpinner(text="Auto",
                                                     values=["Auto", "PCR", "Anneal"],
                                                     pos_hint={"x": self.parent.dropdown_width * self.my_id + self.parent.xpos,
                                                               "top": self.parent.ypos - self.parent.dropdown_height},
                                                     size_hint=(self.parent.dropdown_width, self.parent.dropdown_height),
                                                     # background_color=[0.898, 0.584, 0, 1],
                                                     # background_normal='primer_spinner_up.png',
                                                     # background_down='primer_spinner_down.png',
                                                     background_normal='primer_spin_up.png',
                                                     background_down='primer_spin_up.png',
                                                     option_cls=MySpinnerOption,
                                                     my_id=self.my_id,
                                                     my_row=self.my_row))
            if own_number not in oh_ids and own_number < len(selection) - 1 and own_number != 0:
                if len(selection) > own_number + 1 and not any([True for frag in selection[own_number + 1] if registered_fragments[frag][1] in ["plasmid", "dprom", "dseed", "dscaf"]]):
                    self.parent.add_widget(OhInput(text=f"overhang",
                                                   size_hint=(self.parent.dropdown_width, self.parent.dropdown_height),
                                                   pos_hint={"x":  0.5 * self.parent.dropdown_width + self.parent.dropdown_width * self.my_id + self.parent.xpos,
                                                             "top": self.parent.ypos - (2 * self.parent.dropdown_height) + (self.parent.row_count - self.my_row) * (self.parent.dropdown_height + 0.1)},
                                                   my_id=self.my_id,
                                                   my_row=self.my_row))
            if own_number - 1 not in oh_ids and own_number - 1 != 0:
                if not any([True for frag in selection[own_number - 1] if registered_fragments[frag][1] in ["plasmid", "dprom", "dseed", "dscaf"]]):
                    self.parent.add_widget(OhInput(text=f"overhang",
                                                   size_hint=(self.parent.dropdown_width, self.parent.dropdown_height),
                                                   pos_hint={"x": 0.5 * self.parent.dropdown_width + self.parent.dropdown_width * prev_id + self.parent.xpos,
                                                             "top": self.parent.ypos - (2 * self.parent.dropdown_height) + (self.parent.row_count - prev_row) * (self.parent.dropdown_height + 0.1)},
                                                   my_id=prev_id,
                                                   my_row=prev_row))


########################################################################################################################
# Combination of TextInput and submit Button for adding single fragments to saved_data.csv #######################
########################################################################################################################


class SubmitText(Widget):

    inputs = StringProperty("")

    def __init__(self, **kwargs):
        super(SubmitText, self).__init__(**kwargs)
        self.bind(inputs=lambda x, seq: self.check_fragment_seq())

    def check_fragment_seq(self):
        allowed = ["A", "C", "G", "T"]
        allowed_frags = ["plasmid", "promoter", "scaffold", "seed", "dprom", "dscaf", "dseed"]
        split_inputs = self.inputs.split(",")
        if not all([base in allowed for base in split_inputs[2].upper()]) and len(self.inputs) > 0:
            # print(self.inputs)
            raise ValueError("There are unallowed bases in your input sequence")
        elif split_inputs[1] not in allowed_frags:
            raise ValueError("The fragment type you specified doesn't exist")
        else:
            self.add_fragment()

    def add_fragment(self):
        split_inputs = self.inputs.split(",")
        split_inputs[2] = split_inputs[2].upper()
        print(split_inputs)
        registered_fragments[split_inputs[0]] = split_inputs
        with open('data/saved_data.csv', "a", newline="") as csvfile:
            content_adder = csv.writer(csvfile, delimiter=";")
            content_adder.writerow(split_inputs)


########################################################################################################################
# Predefined SpinnerOption parameters for PrimerSpinner ################################################################
########################################################################################################################


class MySpinnerOption(SpinnerOption):
    def __init__(self, **kwargs):
        super(MySpinnerOption, self).__init__(**kwargs)
        # self.background_color = [0.898, 0.584, 0, 1]
        self.background_down = "primer_spin_up.png"
        self.background_normal = "primer_spin_up.png"



########################################################################################################################
# Wrapper class for the dropdown construct consisting of AddButton, MultiSelectSpinner, OhInput and PrimerSpinner ######
########################################################################################################################


class MultiMulti(FloatLayout):

    xpos = NumericProperty(0.01)

    ypos = NumericProperty(0.99)

    dropdown_width = NumericProperty(0.1)

    dropdown_height = NumericProperty(0.05)

    adder_width = NumericProperty(0.05)

    adder_height = NumericProperty(0.05)

    drop_count = NumericProperty(0)

    row_count = NumericProperty(0)

    run_input_fragments = []

    def __init__(self, **kwargs):
        super(MultiMulti, self).__init__(**kwargs)
        self.dropdown_keys = self.get_values()
        # print(self.dropdown_keys)

    def add_dropdown(self, add_button_id, del_button_id, get_seq_button_id):
        if (self.drop_count + 2) * self.dropdown_width % 1 == 0:
            self.ypos -= self.dropdown_height + 0.1
            self.drop_count = 0
            self.row_count += 1

        plasmid = True if self.drop_count == 0 and self.row_count == 0 else False

        self.update_values(plasmid=plasmid)
        self.add_widget(MultiSelectSpinner(text="",
                                           values=self.dropdown_keys,
                                           pos_hint={"x": self.dropdown_width * self.drop_count + self.xpos,
                                                     "top": self.ypos},
                                           size_hint=(self.dropdown_width, self.dropdown_height),
                                           my_id=self.drop_count,
                                           my_row=self.row_count))
        self.add_widget(PrimerSpinner(text="Auto",
                                      values=["Auto", "PCR", "Anneal"],
                                      pos_hint={"x": self.dropdown_width * self.drop_count + self.xpos,
                                                "top": self.ypos - self.dropdown_height},
                                      size_hint=(self.dropdown_width, self.dropdown_height),
                                      # background_color=[0.898, 0.584, 0, 1],
                                      # background_normal='primer_spinner_up.png',
                                      # background_down='primer_spinner_down.png',
                                      background_normal='primer_spin_up.png',
                                      background_down='primer_spin_up.png',
                                      option_cls=MySpinnerOption,
                                      my_id=self.drop_count,
                                      my_row=self.row_count))

        add_button_id.pos_hint = {"x": self.dropdown_width * (self.drop_count + 1) + self.xpos, "top": self.ypos}
        del_button_id.pos_hint = {"x": self.dropdown_width * (self.drop_count + 1) + self.xpos, "top": self.ypos - 0.05}
        get_seq_button_id.pos_hint = {"x": 0.05 + self.dropdown_width * (self.drop_count + 1) + self.xpos, "top": self.ypos - 0.025}
        self.drop_count += 1

        if self.drop_count == 0 and self.row_count == 0:
            self.check_saved_data()

    def update_values(self, plasmid=False):
        dropdown_keys = []
        with open('data/saved_data.csv') as csvfile:
            save_reader = csv.reader(csvfile, delimiter=";")
            if not plasmid:
                for idx, row in enumerate(save_reader):
                    if idx > 0 and row[1] != "plasmid":
                        line_key = row[0]
                        dropdown_keys.append(line_key)
                self.dropdown_keys = sorted(dropdown_keys)
            else:
                for idx, row in enumerate(save_reader):
                    if idx > 0 and row[1] == "plasmid":
                        line_key = row[0]
                        dropdown_keys.append(line_key)
                self.dropdown_keys = sorted(dropdown_keys)

    def get_values(self):
        dropdown_keys = []
        with open('data/saved_data.csv') as csvfile:
            save_reader = csv.reader(csvfile, delimiter=";")
            for idx, row in enumerate(save_reader):
                if idx > 0:
                    if row[1] != "plasmid":
                        if len(row) == 3:
                            registered_fragments[row[0]] = row
                        elif len(row) == 4:
                            registered_fragments[row[0]] = row[:-1]
                        else:
                            raise ValueError(f"Fragment {row[0]} contains wrong number of elements")
                    else:
                        registered_fragments[row[0]] = row
                    line_key = row[0]
                    dropdown_keys.append(line_key)
        return dropdown_keys

########################################################################################################################
# Button that is used to produce intermediate_input_dummy.csv and start the tool_main.py script while catching errors ##
########################################################################################################################


class RunButton(Button):

    def create_products(self, selected):
        products = list(itertools.product(*selected))
        uniques = itertools.chain.from_iterable(selected)
        for fragment in uniques:
            check_seq = registered_fragments[fragment][2]
            check_seq = check_seq.replace("[", "")
            check_seq = check_seq.replace("]", "")
            if detect_illegal_res_sites(querry=check_seq, motif=res_2s[self.parent.ids["meta"].meta_dict["Restriction Enzyme"]][0]):
                raise ValueError(f"Fragment {fragment} contains recognition site for type IIS enzyme")
        return products

    def start_protocol(self):

        non_unique_error = False
        oh_collection = []

        try:
            dt_dir = self.get_time()
            if not os.path.isdir(f"data/output/{self.parent.ids['meta'].meta_dict['Project Name']}_{dt_dir}"):
                os.mkdir(f"data/output/{self.parent.ids['meta'].meta_dict['Project Name']}_{dt_dir}")
            else:
                raise ValueError("Date-time directory already exists")
            if not os.path.isdir(f"data/interm_input/{self.parent.ids['meta'].meta_dict['Project Name']}_{dt_dir}"):
                os.mkdir(f"data/interm_input/{self.parent.ids['meta'].meta_dict['Project Name']}_{dt_dir}")
            else:
                raise ValueError("Date-time directory already exists")

            coll = open(f"data/output/{self.parent.ids['meta'].meta_dict['Project Name']}_{dt_dir}/{self.parent.ids['meta'].meta_dict['Project Name']}_coll.csv", mode="w", newline="")
            content_adder = csv.writer(coll, delimiter=";")
            content_adder.writerow(["Product No.", "Primer", "Sequence (5' to 3')"])
            coll.close()

            start = 1
            meta_dummy = []

            for meta_par in self.parent.ids["meta"].meta_dict:
                check_meta_inputs(id=meta_par, value=self.parent.ids['meta'].meta_dict[meta_par])
                if meta_par == "Project Name":
                    meta_dummy.append("{}={}".format(translations[meta_par], f"data/output/{self.parent.ids['meta'].meta_dict['Project Name']}_{dt_dir}/{self.parent.ids['meta'].meta_dict[meta_par]}"))
                else:
                    meta_dummy.append("{}={}".format(translations[meta_par], self.parent.ids["meta"].meta_dict[meta_par]))

            products = self.create_products(selected=selection)
            if len(products) > 0 and selection != []:
                for idx, product in enumerate(products):
                    print(f"product = {product}")
                    if len(registered_fragments[product[0]]) == 4 and registered_fragments[product[0]][1] == "plasmid":
                        meta_dummy.append(f"plas_gb={registered_fragments[product[0]][3]}")
                        print(meta_dummy)
                    else:
                        meta_dummy.append("plas_gb=None")
                        print(meta_dummy)
                    with open(f"./data/interm_input/{self.parent.ids['meta'].meta_dict['Project Name']}_{dt_dir}/intermediate_input_dummy_{idx + int(start)}.csv", "w", newline="") as csv_input:
                        csv_writer = csv.writer(csv_input)
                        csv_writer.writerow(meta_dummy)
                        for idx2, key in enumerate(product):
                            # list has to be sliced in order to create a copy of the list from registered fragments
                            # otherwise changes in fragment will be applied to registered fragments
                            fragment = registered_fragments[key][:]
                            fragment.append(custom_primer[idx2])
                            fragment.append(custom_ohs[idx2])
                            csv_writer.writerow(fragment)
                    print("Running tool")
                    prims, ohs = all_in_one_better(in_file="./data/interm_input/{}_{}/intermediate_input_dummy_{}.csv".format(self.parent.ids['meta'].meta_dict['Project Name'], dt_dir,  idx + int(start)), constr_nr=(idx + int(start)))
                    oh_collection.append(ohs)
                    if prims == "non_unique":
                        non_unique_error = True
                if non_unique_error:
                    # self.parent.ids["vis_out"].text = f"[ref=reset]Detected duplicates in overhangs[/ref]"
                    oh_view = OhViewPopup(ohs=oh_collection, prods=products, mat_file="FileS03_T4_18h_25C.xlsx")
                    print(f"SELECTION={selection}")
                    print(f"OH COLLECTION={oh_collection}")
                    oh_view.open()
                    if os.path.isdir(f"data/output/{self.parent.ids['meta'].meta_dict['Project Name']}_{dt_dir}"):
                        shutil.rmtree(f"data/output/{self.parent.ids['meta'].meta_dict['Project Name']}_{dt_dir}")
                    if os.path.isdir(f"data/interm_input/{self.parent.ids['meta'].meta_dict['Project Name']}_{dt_dir}"):
                        shutil.rmtree(f"data/interm_input/{self.parent.ids['meta'].meta_dict['Project Name']}_{dt_dir}")
                    raise ValueError("Detected duplicates in overhangs")
                else:
                    # print("started oh opening")
                    # oh_view = OhViewPopup(ohs=oh_collection, prods=products)
                    # oh_view.open()
                    pass
            else:
                raise ValueError("No products were returned\nCheck if all Fragments are selected")

            print(registered_fragments)
            build_prim_coll(f"data/output/{self.parent.ids['meta'].meta_dict['Project Name']}_{dt_dir}/{self.parent.ids['meta'].meta_dict['Project Name']}")
            # self.get_root_window().add_widget(Popup())

        except Exception as exc:
            if os.path.isdir(f"data/output/{dt_dir}_{self.parent.ids['meta'].meta_dict['Project Name']}"):
                shutil.rmtree(f"data/output/{dt_dir}_{self.parent.ids['meta'].meta_dict['Project Name']}")
            if os.path.isdir(f"data/interm_input/{dt_dir}_{self.parent.ids['meta'].meta_dict['Project Name']}"):
                shutil.rmtree(f"data/interm_input/{dt_dir}_{self.parent.ids['meta'].meta_dict['Project Name']}")
            self.parent.ids["vis_out"].text = f"[ref=reset]{str(exc)}[/ref]"

    def add_linklabel(self, parent):
        present = False
        for child in parent.children:
            if "LinkLabel" in str(child):
                present = True
        if not present:
            my_links = LinkLabel(pos=(parent.width, parent.height))
            parent.add_widget(my_links)

    def get_time(self):
        time = str(datetime.datetime.now())
        time = time.split(".")[0]
        time = time.replace("-", "_")
        time = time.replace(":", "_")
        time = time.replace(" ", "_")
        return time

########################################################################################################################
# Popup that is used to select a fasta file containing fragments which should be added to saved_data.csv #########
########################################################################################################################


class MyPopup(Popup):

    def print_cont(self, file):
        if file[-6:] == ".fasta" or file[-3:] == ".fa":
            content_to_add = read_in_file_fasta(file)
        elif file[-4:] == ".csv":
            content_to_add = read_seq_from_csv(file)
        else:
            raise ValueError("Trying to load from file with incorrect ending")
        c = open("data/saved_data.csv", "a", newline="")
        content_adder = csv.writer(c, delimiter=";")
        # content_adder.writerow([])
        for line in content_to_add:
            if line[0] not in registered_fragments.keys():
                registered_fragments[line[0]] = line
                content_adder.writerow(line)
            else:
                print("Element key {} from input file is already used".format(line[0]))
        c.close()


class OhViewPopup(Popup):
    def __init__(self, ohs, prods, mat_file=False, **kwargs):
        super(OhViewPopup, self).__init__(**kwargs)
        self.ohs = ohs
        self.prods = prods
        self.mat_file = mat_file
        print(self.ohs)
        print(self.prods)
        try:
            doubles = [[oh for oh_pair in construct for oh in oh_pair] for construct in self.ohs]
            # print(doubles)
            # doubles = [[oh_pair[1] for oh_pair in construct] for construct in self.ohs]
            counts = [Counter(double) for double in doubles]

            start_y = 1.45
            set_fid = -1
            app = App.get_running_app()
            if self.mat_file:
                excel_data = pd.read_excel(f'sb8b00333_si_002/Supplemental Data {app.root.ids["meta"].ids["oh_len"].ids["txt"].text}/{self.mat_file}', header=0, index_col=0)
                data = pd.DataFrame(excel_data)
            for mul_idx, mul in enumerate(self.prods):

                oh_fids = {"green": [], "red": [], "yellow": []}
                prod_string = ""
                print(f"OVERHANGS = {list(itertools.chain(*ohs[mul_idx]))}")
                start_x = 0
                if len(ohs[mul_idx]) > 0 and self.mat_file:
                    fids, set_fid = ligation_fidelity(ohs=list(itertools.chain(*ohs[mul_idx])), mat=data)
                    print(f"FIDS = {fids}")
                p_ges = 1
                for frag_idx, frag in enumerate(mul):
                    if self.mat_file:
                        p_oh_pair = fids[self.ohs[mul_idx][frag_idx][1]]
                        c_oh_pair = "f3240c" if int(counts[mul_idx][self.ohs[mul_idx][frag_idx][1]]) > 1 else "50ea23" if (p_oh_pair >= 0.99 and fids[self.ohs[mul_idx][frag_idx][0]] >= 0.99) else "feec08"
                        p_ges *= p_oh_pair

                        if c_oh_pair == "f3240c":
                            oh_fids["red"].append(ohs[mul_idx][frag_idx][1])
                        elif c_oh_pair == "feec08":
                            oh_fids["yellow"].append(ohs[mul_idx][frag_idx][1])
                        else:
                            oh_fids["green"].append(ohs[mul_idx][frag_idx][1])

                    prod_string += f"{frag}  "

                    prod_string += f"[color=#{c_oh_pair}]-{self.ohs[mul_idx][frag_idx][1]}-[/color]  "

                if set_fid != -1:

                    prod_string += f"Fidelity={set_fid}"

                lab_fid = OhViewLabel(
                    text=f"[ref=cp_ohs]{prod_string}[/ref]",
                    pos_hint={"center_y": start_y, "x": start_x},
                    ohs=self.ohs[mul_idx],
                    outline_width=5,
                    counts=counts[mul_idx],
                    oh_fids=oh_fids
                )

                print(prod_string)
                print(oh_fids)

                self.ids["canvas"].add_widget(lab_fid)

                start_y -= 0.05
                # print(f"p_ges={p_ges}")
        except Exception as exc:
            self.parent.ids["vis_out"].text = f"[ref=reset]{str(exc)}[/ref]"

    def update_fids(self, mat):
        file_conns = {"T4 01h 25C": "FileS01_T4_01h_25C.xlsx",
                      "T4 01h 37C": "FileS02_T4_01h_37C.xlsx",
                      "T4 18h 25C": "FileS03_T4_18h_25C.xlsx",
                      "T4 18h 37C": "FileS04_T4_18h_37C.xlsx"}

        app = App.get_running_app()
        print(app.root.ids["meta"].ids["oh_len"].ids["txt"].text)

        excel_data = pd.read_excel(f'sb8b00333_si_002/Supplemental Data {app.root.ids["meta"].ids["oh_len"].ids["txt"].text}/{file_conns[mat]}', header=0, index_col=0)
        data = pd.DataFrame(excel_data)

        for widget in self.ids["canvas"].children:
            if "OhViewLabel" in str(widget):
                oh_count = 0
                oh_fids = {"green": [], "red": [], "yellow": []}
                widget_content = widget.text
                widget_content = widget_content.split("  ")
                fid = widget_content[-1].split("=")
                frag_ohs = widget_content[:-1]
                fids, set_fid = ligation_fidelity(ohs=list(itertools.chain(*widget.ohs)), mat=data)
                fid[1] = str(set_fid)
                fid = "=".join(fid)
                for info_id, information in enumerate(frag_ohs):
                    if "color=" in information:
                        c_oh_pair = "f3240c" if int(widget.counts[widget.ohs[oh_count][1]]) > 1 else "50ea23" if (fids[widget.ohs[oh_count][1]] >= 0.99 and fids[widget.ohs[oh_count][0]] >= 0.99) else "feec08"
                        frag_ohs[info_id] = f"{frag_ohs[info_id][:8]}{c_oh_pair}{frag_ohs[info_id][14:]}"
                        oh_count += 1
                        if c_oh_pair == "f3240c":
                            oh_fids["red"].append(frag_ohs[info_id][16:16+int(app.root.ids["meta"].ids["oh_len"].ids["txt"].text)])
                        elif c_oh_pair == "feec08":
                            oh_fids["yellow"].append(frag_ohs[info_id][16:16+int(app.root.ids["meta"].ids["oh_len"].ids["txt"].text)])
                        else:
                            oh_fids["green"].append(frag_ohs[info_id][16:16+int(app.root.ids["meta"].ids["oh_len"].ids["txt"].text)])
                widget.text = "  ".join(frag_ohs + [fid])
                widget.oh_fids = oh_fids


class FormatPopup(Popup):
    pass


class SliderPopup(Popup):
    pass


class FormatButton(Button):
    def open_format(self):
        format_view = FormatPopup()
        format_view.open()


class OhViewLabel(Label):
    position = NumericProperty(0)
    counts = {}
    ohs = []
    oh_fids = {}

    def __init__(self, oh_fids, counts=None, ohs=None, **kwargs):
        super(OhViewLabel, self).__init__(**kwargs)
        self.ohs = ohs
        self.counts = counts
        self.oh_fids = oh_fids

    def cp_ohs(self):
        col = np.random.random(4)
        col[3] = 1
        print(col)
        self.outline_color = col
        sel_complexity = self.parent.parent.parent.parent.parent.ids["selection_slider"].value
        if sel_complexity == 0:
            pyperclip3.copy(f"{self.ohs[0][1]},{self.ohs[-1][1]}")
            print(f"{self.ohs[0][1]},{self.ohs[-1][1]}")
        elif sel_complexity == 1:
            pyperclip3.copy(f"{','.join(self.oh_fids['green'])}")
            print(f"{','.join(self.oh_fids['green'])}")
        elif sel_complexity == 2:
            pyperclip3.copy(f"{','.join(self.oh_fids['green'])},{','.join(self.oh_fids['yellow'])}")
            print(f"{','.join(self.oh_fids['green'])},{','.join(self.oh_fids['yellow'])}")
        else:
            raise Exception("Slider state not allowed")


class GetSeq(Button):
    def open_get_seq(self):
        get_seq_view = GetSeqPopup()
        get_seq_view.open()


class GetSeqPopup(Popup):
    def __init__(self, **kwargs):
        super(GetSeqPopup, self).__init__(**kwargs)
        fragments = list(itertools.chain(*selection))
        fragments = np.unique(fragments)
        fragments.sort()
        print(fragments)
        start_y_center = 1.455
        start_y_top = 0.96
        for fragment in fragments:
            # if registered_fragments[fragment][1] != "plasmid":
            cp_button = CopyButton(text="copy",
                                   pos_hint={"x": 0.025, "center_y": start_y_top},
                                   size_hint=(0.075, 0.02),
                                   sequence=registered_fragments[fragment][2])
            frag_seq = SeqLabel(text=f"{fragment} - {registered_fragments[fragment][2]}",
                                pos_hint={"x": 0.04, "center_y": start_y_center})
            self.ids["canvas"].add_widget(cp_button)
            self.ids["canvas"].add_widget(frag_seq)
            start_y_center -= 0.025
            start_y_top -= 0.025


class CopyButton(Button):
    sequence = StringProperty("")

    def copy_seq(self):
        print(self.sequence)
        pyperclip3.copy(self.sequence)


class SeqLabel(Label):
    pass


class FancySlider(Slider):
    def __init__(self, **kwargs):
        super(FancySlider, self).__init__(**kwargs)
        self.bind(value=self.on_slide)

    def on_slide(self, instance, value):
        if self.value == 0:
            self.cursor_image = "./slider_img_white.png"
        elif self.value == 1:
            self.cursor_image = "./slider_img_red.png"
        elif self.value == 2:
            self.cursor_image = "./slider_img_red_yellow.png"


class SliderButton(Button):
    def open_format(self):
        slider_view = SliderPopup()
        slider_view.open()

########################################################################################################################
# App Wrapper ##########################################################################################################
########################################################################################################################


class FinalGuiApp(App):

    def build(self):
        self.title = 'G-GArden'


if __name__ == "__main__":
    FinalGuiApp().run()
