import subprocess

# import inspect
# import sys

import papermill as pm

# sys.path.append(f"{Path(inspect.getfile(inspect.currentframe())).absolute().parent.parent}")
from General.type_hints import StrOrPath


def execute_notebook(
    template_notebook: StrOrPath,
    execution_out_dir: StrOrPath,
    executed_notebook_base_name: str,
    **parameters,
):

    full_executed_notebook_ipynb_path = (
        f"{execution_out_dir}/{executed_notebook_base_name}.ipynb"
    )
    clean_executed_notebook_html_path = (
        f"{execution_out_dir}/{executed_notebook_base_name}.no-input.no-prompt.html"
    )

    print(f"{template_notebook = }")
    print(f"{execution_out_dir = }")
    print(f"{executed_notebook_base_name = }")
    print(f"{parameters = }")
    print(f"{full_executed_notebook_ipynb_path = }")
    print(f"{clean_executed_notebook_html_path = }")
    # return   # todo delete

    pm.execute_notebook(
        template_notebook, full_executed_notebook_ipynb_path, parameters=parameters
    )

    full_convert_cmd = (
        f"jupyter nbconvert --to html {full_executed_notebook_ipynb_path} "
        f"--TagRemovePreprocessor.remove_input_tags item"
    )
    subprocess.run(full_convert_cmd, shell=True)

    clean_convert_cmd = (
        f"jupyter nbconvert --to html --stdout --no-input --no-prompt "
        f"{full_executed_notebook_ipynb_path} > {clean_executed_notebook_html_path}"
    )
    subprocess.run(clean_convert_cmd, shell=True)
