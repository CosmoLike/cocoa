#! /usr/bin/env python
# encoding: utf-8
# WARNING! Do not edit! https://waf.io/book/index.html#_obtaining_the_waf_file

import re
from waflib import Task,Logs
from waflib.TaskGen import extension
cy_api_pat=re.compile(r'\s*?cdef\s*?(public|api)\w*')
re_cyt=re.compile(r"""
	^\s*                           # must begin with some whitespace characters
	(?:from\s+(\w+)(?:\.\w+)*\s+)? # optionally match "from foo(.baz)" and capture foo
	c?import\s(\w+|[*])            # require "import bar" and capture bar
	""",re.M|re.VERBOSE)
@extension('.pyx')
def add_cython_file(self,node):
	ext='.c'
	if'cxx'in self.features:
		self.env.append_unique('CYTHONFLAGS','--cplus')
		ext='.cc'
	for x in getattr(self,'cython_includes',[]):
		d=self.path.find_dir(x)
		if d:
			self.env.append_unique('CYTHONFLAGS','-I%s'%d.abspath())
	tsk=self.create_task('cython',node,node.change_ext(ext))
	self.source+=tsk.outputs
class cython(Task.Task):
	run_str='${CYTHON} ${CYTHONFLAGS} -o ${TGT[0].abspath()} ${SRC}'
	color='GREEN'
	vars=['INCLUDES']
	ext_out=['.h']
	def runnable_status(self):
		ret=super(cython,self).runnable_status()
		if ret==Task.ASK_LATER:
			return ret
		for x in self.generator.bld.raw_deps[self.uid()]:
			if x.startswith('header:'):
				self.outputs.append(self.inputs[0].parent.find_or_declare(x.replace('header:','')))
		return super(cython,self).runnable_status()
	def post_run(self):
		for x in self.outputs:
			if x.name.endswith('.h'):
				if not x.exists():
					if Logs.verbose:
						Logs.warn('Expected %r',x.abspath())
					x.write('')
		return Task.Task.post_run(self)
	def scan(self):
		node=self.inputs[0]
		txt=node.read()
		mods=set()
		for m in re_cyt.finditer(txt):
			if m.group(1):
				mods.add(m.group(1))
			else:
				mods.add(m.group(2))
		Logs.debug('cython: mods %r',mods)
		incs=getattr(self.generator,'cython_includes',[])
		incs=[self.generator.path.find_dir(x)for x in incs]
		incs.append(node.parent)
		found=[]
		missing=[]
		for x in sorted(mods):
			for y in incs:
				k=y.find_resource(x+'.pxd')
				if k:
					found.append(k)
					break
			else:
				missing.append(x)
		implicit=node.parent.find_resource(node.name[:-3]+'pxd')
		if implicit:
			found.append(implicit)
		Logs.debug('cython: found %r',found)
		has_api=False
		has_public=False
		for l in txt.splitlines():
			if cy_api_pat.match(l):
				if' api 'in l:
					has_api=True
				if' public 'in l:
					has_public=True
		name=node.name.replace('.pyx','')
		if has_api:
			missing.append('header:%s_api.h'%name)
		if has_public:
			missing.append('header:%s.h'%name)
		return(found,missing)
def options(ctx):
	ctx.add_option('--cython-flags',action='store',default='',help='space separated list of flags to pass to cython')
def configure(ctx):
	if not ctx.env.CC and not ctx.env.CXX:
		ctx.fatal('Load a C/C++ compiler first')
	if not ctx.env.PYTHON:
		ctx.fatal('Load the python tool first!')
	ctx.find_program('cython',var='CYTHON')
	if hasattr(ctx.options,'cython_flags'):
		ctx.env.CYTHONFLAGS=ctx.options.cython_flags
