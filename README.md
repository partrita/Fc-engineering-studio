# Fc Engineering Studio Pro

Fc_Hinge–CH2–CH3 영역을 파라미터화해 charge, glycosylation, FcγR binding hot-spot을 기반으로 설계하는 TUI(Text User Interface) 기반 도구입니다.

이 도구는 **EU Numbering** 체계를 기준으로 인간 IgG Fc 변이 서열을 신속하게 생성하고 검증할 수 있도록 도와줍니다.

## 주요 기능

- **Isotype & Allotype 선택**: IgG1, IgG2, IgG4의 다양한 Allotype(WT, Trastuzumab, Rituximab 등) 베이스 서열 지원
- **Preset Mutations (Common Mutations)**: 널리 알려진 주요 변이들을 체크박스로 간편하게 적용
  - LALA (Effector Silencing)
  - YTE/LS (Half-life extension)
  - Knob-into-Hole (Bispecific design)
  - S228P (IgG4 Hinge stabilization) 등
- **Custom Mutants 입력**: 사용자가 원하는 임의의 변이를 `S239D/I332E`와 같은 형식으로 직접 입력 가능
- **FASTA 생성 및 복사**: 생성된 서열을 즉시 FASTA 형식으로 확인하고 클립보드에 복사

## 설치 및 실행

이 프로젝트는 [uv](https://github.com/astral-sh/uv)를 사용하여 의존성을 관리합니다.

### 실행 방법
```bash
uv run fc-engineer
```

### 개발 환경 설정
```bash
# 의존성 설치
uv sync

# 테스트 실행
uv run pytest
```

## 사용법 및 단축키

1. 왼쪽 패널에서 **Isotype**과 **Allotype**을 선택합니다.
2. 적용하고자 하는 **Common Mutations**를 체크합니다.
3. 추가적인 변이가 있다면 **Custom Mutants** 창에 입력합니다 (슬래시`/` 또는 콤마`,`로 구분).
4. **Generate** 버튼을 누르거나 `Enter` 키를 입력하여 서열을 생성합니다.

### 단축키 (Shortcuts)
- `Enter`: 서열 생성 (Generate)
- `Ctrl + Y`: 결과 서열 클립보드 복사 (Copy)
- `Ctrl + C`: 프로그램 종료 (Quit)

## 프로젝트 구조

- `src/main.py`: TUI 및 핵심 로직 구현
- `src/sequences.yaml`: Isotype별 베이스 서열 데이터
- `src/mutants.yaml`: 주요 변이 프리셋 데이터
- `pyproject.toml`: 프로젝트 설정 및 의존성 정의

## 데이터 커스터마이징

새로운 베이스 서열이나 변이 프리셋을 추가하려면 `src/sequences.yaml` 또는 `src/mutants.yaml` 파일을 수정하십시오. 프로그램 실행 시 해당 파일들을 자동으로 로드합니다.
