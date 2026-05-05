# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 09:00:23 2026

@author: petar.hristov
"""
import os
import ast

REPLACEABLE = ['min', 'max', 'abs', 'sqrt', 'log', 'exp',
               'sin', 'cos', 'tan', 'asin', 'acos', 'atan']

class NpToIvalTransformer(ast.NodeTransformer):
    def __init__(self, whitelist, np_alias="np", new_alias="ival"):
        self.whitelist = set(whitelist)
        self.np_alias = np_alias
        self.new_alias = new_alias
        self.found_numpy_import = False
        self.has_ival_import = False

    def visit_Import(self, node):
        for alias in node.names:
            if alias.name == "numpy" and (alias.asname or "numpy") == self.np_alias:
                self.found_numpy_import = True
            if alias.name == self.new_alias:
                self.has_ival_import = True
        return node

    def visit_ImportFrom(self, node):
        # Optional: detect "from ival import ..." too
        if node.module == self.new_alias:
            self.has_ival_import = True
        return node

    def visit_Call(self, node):
        self.generic_visit(node)

        # Match: np.<func>(...)
        if isinstance(node.func, ast.Attribute):
            attr = node.func

            if (
                isinstance(attr.value, ast.Name)
                and attr.value.id == self.np_alias
                and attr.attr in self.whitelist
            ):
                attr.value.id = self.new_alias

        return node

    def visit_Module(self, node):
        # First process everything
        self.generic_visit(node)

        # If numpy import exists and ival not yet imported → add it
        if self.found_numpy_import and not self.has_ival_import:
            ival_import = ast.Import(names=[ast.alias(name='Interval', asname=self.new_alias)])

            # Insert after last import
            insert_idx = 0
            for i, stmt in enumerate(node.body):
                if isinstance(stmt, (ast.Import, ast.ImportFrom)):
                    insert_idx = i + 1

            node.body.insert(insert_idx, ival_import)

        return node
    
def replace_np_with_ival(code: str, whitelist=REPLACEABLE):
    tree = ast.parse(code)
    transformer = NpToIvalTransformer(whitelist)
    tree = transformer.visit(tree)
    ast.fix_missing_locations(tree)
    return ast.unparse(tree)

def translate(file_name_read, file_name_write=None) -> None:
    #Handle paths
    path_file_read = os.path.dirname(file_name_read)
    file_name_read = os.path.basename(file_name_read)
    
    #Read the file containing the original function
    # imports = 'import Interval as ival\n\n'
    
    with open(file_name_read, 'r') as file:
        # line = file.readline()
        # new_def = imports + replace_np_with_ival(file.read())
        new_def = replace_np_with_ival(file.read())
    
    if not file_name_write:
        split = file_name_read.split('.')
        name, ext = split[:-1], split[-1] #Handle the extension correctly
        file_name_write = ''.join(name) + '_ival.' + ext
        file_name_write = os.path.join(path_file_read, file_name_write)
    
    if os.path.dirname(file_name_read) == '': #New file name has no directory
        file_name_write = os.path.join(path_file_read, file_name_write)
        
    with open(file_name_write, 'w') as file:
        file.write(new_def)

