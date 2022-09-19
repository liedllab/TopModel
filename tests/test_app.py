import pytest
from click.testing import CliRunner
from topmodel import app as application
from topmodel.util.parser import get_structure


file = './data/modelled/output.pdb'
struc = get_structure(file)


def test_output_failure():
    app = application.App(True, True, True)
    with pytest.raises(ValueError):
        app.output_to_terminal()


def test_score_method_failure():
    app = application.App(True, False, False)
    with pytest.raises(ValueError):
        app.compute_score()


def test_score_attr_failure():
    app = application.App(True, False, False)
    with pytest.raises(ValueError):
        app.score


def test_output_exists(capsys):
    app = application.App(True, False, False)
    app.process_structure(struc)
    app.output_to_terminal()
    assert len(capsys.readouterr().out) > 0


@pytest.mark.parametrize('file,exit_code',
                         [[file, 0],
                          ['./data/alanine_dipeptide/1.pdb', 0],
                          ['x,asd', 1],
                          ['./data/fileformats/7SG5_model.pdb', 0]],  # now works
                         )
def test_app_main_correct_imput(file, exit_code):
    runner = CliRunner()
    result = runner.invoke(application.main,
                           args=f"{file} --score",
                           )
    assert len(result.output) > 0
    assert result.exit_code == exit_code
