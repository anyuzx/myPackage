# this module is to provide a wrapper to generate various lammps data file
import sys
import importlib

__all__ = ['data']

class data:
	def __init__(self):
		self.model_dic = {'SPC/E ice 1h':'ice1h','lennard-jones liquid':'lj',\
						  'SPC/E water':'spce','TIP5P water':'tip5p',\
						  'self-avoiding chain':'lattice_chain',\
						  'lennard-jones polymer':'lj_chain','lennard-jones WLC polymer':'lj_wlc',\
						  'Rosette Loop Model polymer':'RSM_chain','Linear Loop Model polymer':'LLM_chain',\
						  'Cross Loop Model polymer':'CLM_chain','Fractal Loop Model polymer':'FLM_chain'}
		
	def available_models(self):
		string = ''
		for key,value in self.model_dic.iteritems():
			string += (key+':').ljust(50)+value+'\n'
		sys.stdout.write(string)

	def writetofile(self,modelname,foutname,*argv,**kwargs):
		self.foutname = foutname

		if modelname not in self.model_dic.values():
			raise Exception('ERROR: Cannot load model named '+modelname+'. Use the available_models class to check all avaialble models.')
		else:
			self.modelname = modelname

		try:
			name = importlib.import_module('.models.'+self.modelname+'_datafile','myPackage.lammps')
		except ImportError:
			raise

		name.writefile(self.foutname,*argv,**kwargs)


