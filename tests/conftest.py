import pytest
from jec4prompt.main import ProcessingState

@pytest.fixture
def mock_state():
    class MockState(ProcessingState):
        def __init__(self):
            self.logger = None
            
    return MockState()