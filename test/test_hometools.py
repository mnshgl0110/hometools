def test_import():
    from hometools.hometools import readfasta
    assert True

def test_cli():
    from subprocess import Popen, PIPE
    p=Popen('hometools -h'.split(), stdout=PIPE, stderr=PIPE)
    out=p.communicate()
    assert 'error' not in out[1].decode()