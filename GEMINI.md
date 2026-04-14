# Fc Engineering Studio Pro - Project Context

## Project Overview
Fc Engineering Studio Pro는 인간 IgG Fc(Hinge-CH2-CH3) 영역의 변이 서열을 설계하고 분석하는 Python 기반 TUI(Text User Interface) 애플리케이션입니다. **EU Numbering** 체계를 기준으로 단백질 서열의 특정 잔기(Residue) 위치를 계산하고, 변이(Mutation)를 적용하여 FASTA 형식의 최종 서열을 생성합니다.

### 핵심 기술 스택
- **Language**: Python >= 3.14
- **TUI Framework**: [Textual](https://textual.textualize.io/) (Rich 기반의 터미널 인터페이스)
- **Dependency Management**: [uv](https://github.com/astral-sh/uv)
- **Data Persistence**: YAML (`pyyaml`)
- **Clipboard Management**: `pyperclip`
- **Testing**: `pytest`

### 아키텍처
- `src/main.py`: 애플리케이션의 핵심 로직과 TUI 클래스(`MutantApp`)를 포함합니다.
- `src/sequences.yaml`: 각 IgG Isotype(IgG1, IgG2, IgG4) 및 Allotype(WT, Trastuzumab 등)별 기준 서열 데이터를 관리합니다.
- `src/mutants.yaml`: LALA, YTE, LS, Knob-into-Hole 등 업계에서 널리 사용되는 주요 변이 프리셋을 정의합니다.
- `tests/`: 변이 적용 로직 및 번호 체계 계산 로직에 대한 단위 테스트를 포함합니다.

## Building and Running

### 실행 명령어
```bash
# 애플리케이션 실행
uv run fc-engineer

# 또는 진입점 직접 실행
uv run python src/main.py
```

### 개발 및 테스트
```bash
# 의존성 동기화
uv sync

# 테스트 실행
uv run pytest
```

## Development Conventions

### Coding Style
- **Python Typing**: `typing` 모듈을 사용하여 적극적으로 타입을 명시합니다 (`Dict`, `List`, `Optional`, `Tuple` 등).
- **Docstrings**: 클래스 및 주요 함수에 대한 설명을 포함하여 코드 가독성을 유지합니다.
- **Entry Points**: `[project.scripts]`를 통해 터미널 명령어를 제공하며, `main.py` 내부에 `main()` 함수를 두어 실행을 래핑합니다.

### Testing Practices
- `pytest`를 사용하여 로직의 정확성을 검증합니다.
- 특히 EU Numbering에 따른 인덱스 오프셋 계산(`get_residue_index`)과 변이 문자열 파싱(`parse_mutation`)에 대한 엣지 케이스 테스트를 중시합니다.

### Data Management
- 새로운 서열이나 변이 프리셋은 `src/*.yaml` 파일에 선언적으로 추가합니다.
- `pyproject.toml`의 `[tool.setuptools.package-data]` 설정을 통해 패키징 시 YAML 파일이 누락되지 않도록 관리합니다.

## User Interface Guidelines
- **Shortcuts**: 사용자의 편의를 위해 `Enter`(생성), `Ctrl+Y`(클립보드 복사), `Ctrl+C`(종료) 등의 단축키를 적극 활용합니다.
- **Notification**: 변이 적용 중 오류가 발생할 경우 TUI 하단의 `Log` 위젯 또는 `notify` 알림을 통해 사용자에게 즉각 피드백을 제공합니다.
