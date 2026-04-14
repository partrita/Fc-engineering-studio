import pytest
from main import MutantApp

@pytest.mark.asyncio
async def test_app_compose():
    app = MutantApp()
    async with app.run_test() as pilot:
        assert app.title == "Fc Engineering Studio Pro"
        # If it reaches here, compose() didn't crash
