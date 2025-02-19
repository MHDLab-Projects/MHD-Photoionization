source .env

#TODO: setup.py needs to be replaced with pyproject.toml
python setup.py develop # For pi_paper_utils package

git submodule update --init

pushd mhdlab
pip install -e .
popd

pushd mhdlab/mhdlab/fileio
git submodule update --init
popd

if [ "$MUNGE_DATA" = true ] ; then
    # Only post process if we need to munge the data, not currently releasing post processing with paper

    pushd PostProcessor
    # This dummy directory is needed to get test discovery to work
    # the 'munged' data folder is really the test data for the post processed data in this repository. 
    # Post processor tests should not run. 
    mkdir -p tests/test_data

    pip install -e .
    popd
fi