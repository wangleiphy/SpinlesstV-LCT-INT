# Copyright 2013 ETHZ.ch Lukas Gamper <lukas.gamper@usystems.ch>
import numpy as np
import tables as tbl
from collections import Set

class nodeList(Set):

    def __init__(self, nodes):
        self.elements = [node._v_name for node in nodes]

    def __str__(self):
        return str(self.elements)

    def __iter__(self):
        return iter(self.elements)

    def __contains__(self, value):
         return value in self.elements

    def __len__(self):
        return len(self.elements)

class archive:

    def __init__(self, filename, mode = 'r'):
        self.file = tbl.openFile(filename, mode)

    @property
    def isopen(self):
        return self.file.isopen == 1

    def close(self):
        self.file.close()

    def __enter__(self):
        if not self.isopen:
            raise tbl.exceptions.ClosedFileError('HDF5 file is closed')
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __setitem__(self, path, value):
        segments = self.__splitPath(path)
        group = self.file.root
        for segment in segments[:-1]:
            try:
                group = self.file.getNode(group, name=segment, classname='Group')
            except tbl.exceptions.NoSuchNodeError:
                group = self.file.createGroup(group, segment)

            try:
                self.file.createArray(group, segments[-1], value)
            except tbl.exceptions.NodeError:
                self.file.removeNode(group, name=segments[-1], recursive=True)
                self.file.createArray(group, segments[-1], value)

    def __getitem__(self, path):
        if path == '/':
            return nodeList(self.file.listNodes(self.file.root))
        segments = self.__splitPath(path)
        group = self.file.root
        for segment in segments[:-1]:
            group = self.file.getNode(group, name=segment, classname='Group')
        try:
            return self.file.getNode(group, name=segments[-1], classname='Array').read()
        except tbl.exceptions.NoSuchNodeError:
            return nodeList(self.file.listNodes(self.file.getNode(group, name=segments[-1], classname='Group')))

    def __splitPath(self, path):
        segments = path.split('/')
        if len(segments[0]) > 0 or min([len(segment) for segment in segments[1:]]) == 0:
            raise tbl.exceptions.NodeError('Invalid path ' + path)
        return segments[1:]
