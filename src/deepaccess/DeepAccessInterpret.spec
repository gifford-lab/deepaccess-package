# -*- mode: python ; coding: utf-8 -*-

block_cipher = None


a = Analysis(['DeepAccessInterpret.py'],
             pathex=['/data/gl/g1/jhammelm/projects/DeepAccessSquared'],
             binaries=[('/data/gl/g1/jhammelm/env/miniconda3/envs/deepaccessaccess/lib/libmkl_core.so','.'),('/data/gl/g1/jhammelm/env/miniconda3/envs/deepaccessaccess/lib/libmkl_intel_thread.so','.'),
                       ('/data/gl/g1/jhammelm/env/miniconda3/envs/deepaccessaccess/lib/libmkl_intel_lp64.so','.'),('/data/gl/g1/jhammelm/env/miniconda3/lib/libpython3.7m.so.1.0','.'),
                       ('/data/gl/g1/jhammelm/env/miniconda3/envs/deepaccessaccess/lib/libmkl_mc3.so','.'),('/data/gl/g1/jhammelm/env/miniconda3/envs/deepaccessaccess/lib/libmkl_def.so','.'),
                       ('/data/gl/g1/jhammelm/env/miniconda3/lib/libpython3.7m.so','.'),
                       ('/data/gl/g1/jhammelm/env/miniconda3/lib/libcblas.so.3','.'),
                       ('/data/gl/g1/jhammelm/env/miniconda3/envs/deepaccessaccess/lib/python3.7/site-packages/pysam/libctabixproxies.cpython-37m-x86_64-linux-gnu.so','.')],
             datas=[('/afs/csail.mit.edu/u/j/jhammelm/.keras/keras.json' ,'.')],
             hiddenimports=['numpy.random.common','scipy.spatial.transform._rotation_groups','_sysconfigdata_m_linux_x86_64-linux-gnu'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='DeepAccessInterpret',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True )

