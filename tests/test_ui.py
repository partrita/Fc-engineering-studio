import pytest
from fc_engineer.app import MutantApp, ResultScreen

@pytest.mark.asyncio
async def test_app_compose():
    app = MutantApp()
    async with app.run_test() as pilot:
        assert app.title == "Fc Engineering Studio"
        # If it reaches here, compose() didn't crash

@pytest.mark.asyncio
async def test_result_screen_renders_diff():
    app = MutantApp()
    async with app.run_test() as pilot:
        app.selected_isotype = "igg1"
        app.selected_allotype = "WT(P01857-1)"
        app.all_mutants = "L234A/L235A"
        await app.push_screen(ResultScreen())
        await pilot.pause()
        # FASTA 생성 및 diff 뷰 렌더링이 크래시 없이 완료되어야 함
        assert app.last_fasta.startswith(">IGG1_Wt(p01857-1)_L234A_L235A")
        mut_seq = app.last_fasta.split("\n")[1]
        assert mut_seq[116:118] == "AA"  # L234A/L235A 적용 확인 (idx 116=EU234)
