import re
from collections import Counter
import random

random.seed ()

TT_table = {"TTT":"F","TCT":"S","TAT":"Y","TGT":"C","TTC":"F","TCC":"S","TAC":"Y","TGC":"C","TTA":"L","TCA":"S","TAA":"*","TGA":"*","TTG":"L","TCG":"S","TAG":"*","TGG":"W","CTT":"L","CCT":"P","CAT":"H","CGT":"R","CTC":"L","CCC":"P","CAC":"H","CGC":"R","CTA":"L","CCA":"P","CAA":"Q","CGA":"R","CTG":"L","CCG":"P","CAG":"Q","CGG":"R","ATT":"I","ACT":"T","AAT":"N","AGT":"S","ATC":"I","ACC":"T","AAC":"N","AGC":"S","ATA":"I","ACA":"T","AAA":"K","AGA":"R","ATG":"M","ACG":"T","AAG":"K","AGG":"R","GTT":"V","GCT":"A","GAT":"D","GGT":"G","GTC":"V","GCC":"A","GAC":"D","GGC":"G","GTA":"V","GCA":"A","GAA":"E","GGA":"G","GTG":"V","GCG":"A","GAG":"E","GGG":"G","---":"-","NNN":"?"}
nucs = set (['A','C','G','T'])

def TT (codon):
    if codon in TT_table:
        return TT_table[codon]
    return '?'

def newick_parser(nwk_str, bootstrap_values, track_tags=None, optional_starting_tags=None):
  if optional_starting_tags is None:
    optional_starting_tags = {}
    
  clade_stack = []
  automaton_state = 0
  current_node_name = ""
  current_node_attribute = ""
  current_node_annotation = ""
  quote_delimiter = None
  name_quotes = {
    "'": 1,
    '"': 1
  }
  
  def add_new_tree_level():
    new_level = {
    "name": None
    };
    the_parent = clade_stack[len(clade_stack) - 1]
    if (not "children" in the_parent):
      the_parent["children"] = [];
    
    clade_stack.append (new_level);
    the_parent["children"].append(clade_stack[len(clade_stack) - 1]);
    clade_stack[len(clade_stack)-1]["original_child_order"] = len(the_parent["children"])
  

  def finish_node_definition():
    nonlocal current_node_name
    nonlocal current_node_annotation
    nonlocal current_node_attribute
    
    this_node = clade_stack.pop()
    if (bootstrap_values and "children" in this_node):
      this_node["bootstrap_values"] = current_node_name
    else:
      this_node["name"] = current_node_name
    
    this_node["attribute"] = current_node_attribute
    this_node["annotation"] = current_node_annotation
    
    try:
    
      if not 'children' in this_node:
        node_tag = "background"
      for k, v in optional_starting_tags.items():
        if this_node["name"].find (k) >= 0:
          node_tag = v
          break
      else:
        '''
        counts = {}
        node_tag = ""
        for n in this_node['children']:
          counts[n["tag"]] = 1 + (counts[n["tag"]] if n["tag"] in counts  else 0)
        if len (counts) == 1:
          node_tag = list (counts.keys())[0]
        '''
        node_tag = "test"
    
      this_node["tag"] = node_tag
    except Exception as e:
      print ("Exception ", e)
    
    # if track_tags is not None:
    #   track_tags[this_node["name"]] = [this_node["tag"], 'children' in this_node]
    if track_tags is not None:
      track_tags[this_node["name"]] = this_node["tag"]
     
    current_node_name = ""
    current_node_attribute = ""
    current_node_annotation = ""
  

  def generate_error(location):
    return {
      'json': None,
      'error':
        "Unexpected '" +
        nwk_str[location] +
        "' in '" +
        nwk_str[location - 20 : location + 1] +
        "[ERROR HERE]" +
        nwk_str[location + 1 : location + 20] +
        "'"
    }


  tree_json = {
    "name" : "root"
  }
  
  clade_stack.append(tree_json);

  space = re.compile("\s")

  for char_index in range (len(nwk_str)):
    try:
      current_char = nwk_str[char_index]
      if automaton_state == 0:
        #look for the first opening parenthesis
        if (current_char == "("):
          add_new_tree_level()
          automaton_state = 1
      elif automaton_state == 1 or automaton_state == 3:
        #case 1: // name
        #case 3: { // branch length
        #reading name
        if (current_char == ":"):
          automaton_state = 3;
        elif current_char == "," or current_char == ")":
          try:
            finish_node_definition()
            automaton_state = 1
            if (current_char == ","):
              add_new_tree_level()
          except Exception as e:
            return generate_error(char_index)
          
        elif (current_char == "("):
          if len(current_node_name) > 0:
            return generate_error(char_index);
          else:
            add_new_tree_level()
          
        elif (current_char in name_quotes):
          if automaton_state == 1 and len(current_node_name) == 0 and len (current_node_attribute) == 0 and len (current_node_annotation) == 0:
            automaton_state = 2
            quote_delimiter = current_char
            continue
          return generate_error(char_index)
        else:
          if (current_char == "{"):
            if len (current_node_annotation):
              return generate_error(char_index)
            else:
              automaton_state = 4
          else:
            if (automaton_state == 3):
              current_node_attribute += current_char;
            else:
              if (space.search(current_char)):
                continue;
              if (current_char == ";"):
                char_index = len(nwk_str)
                break
            current_node_name += current_char;
      elif automaton_state == 2: 
        # inside a quoted expression
        if (current_char == quote_delimiter):
          if (char_index < len (nwk_str - 1)):
            if (nwk_str[char_index + 1] == quote_delimiter):
              char_index+=1
              current_node_name += quote_delimiter;
              continue;

          quote_delimiter = 0
          automaton_state = 1
          continue
        else:
          current_node_name += current_char;
      elif automaton_state == 4:
        ##inside a comment / attribute
        if (current_char == "}"):
          automaton_state = 3
        else:
          if (current_char == "{"):
            return generate_error(char_index);
          current_node_annotation += current_char;
    except Exception as e:
      return generate_error(char_index);

  if (len (clade_stack) != 1):
    return generate_error(len (nwk_str) - 1);

  if (len (current_node_name)):
    tree_json['name'] = current_node_name;

  return {
    'json': tree_json,
    'error': None
  }


def traverse_tree (node, parent, labels, labeler, composition, subs, leaf_label):
    tag = labeler[node["name"]] if parent else None

    if node["name"] in labels:
        node["label"] = labels[node["name"]]  # what codon is present at the current node?
        if parent:
            diff = 0
            for i, c in enumerate (node["label"]):
                if c in nucs:
                    if parent["label"][i] != c and parent["label"][i] in nucs:
                        diff += 1  # count number of substitutions between parent node and child node
                        
            if diff > 0:
                if not tag in subs:
                    subs [tag] = Counter()
                    
                pt = TT(parent["label"])
                nt = TT(node["label"])
                if pt < nt:
                    sub = pt + ":" + nt
                else:
                    sub = nt + ":" + pt

                subs[tag][sub] += 1
                
    else:
        node["label"] = parent["label"]


    if not "children" in node:
      tag = leaf_label
      if not tag in composition:
        composition[tag] = Counter()
      composition[tag][TT(node["label"])] += 1  # counting how many branches are in the focal clade and how many are in other clades?
        
        
    if "children" in node:
        for c in node["children"]:
           traverse_tree (c, node, labels, labeler, composition, subs, leaf_label) 
    