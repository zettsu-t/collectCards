#include <iostream>
#include <random>
#include <vector>

// 演算精度
enum class Precision {
    PRECISION_MIN,
    PRECISION_BUILTIN_FLOAT = PRECISION_MIN,  // C++ float
    PRECISION_BUILTIN_DOUBLE,                 // C++ double
    PRECISION_VARIABLE_PRECISION_FLOAT,  // 任意精度浮動小数だがコンパイル時に精度を固定する
    PRECISION_MAX_PLUS_ONE,
};

using RandomNumberGen  = std::mt19937;
using RandomNumberSeed = RandomNumberGen::result_type;
using CardNumber       = unsigned int;
using CardNumberArray  = std::vector<CardNumber>;
using TrialCount       = uint32_t;
using DrawCount        = size_t;
using ReportedCount = double;
using ReportedPrice = double;

static_assert(sizeof(TrialCount) <= sizeof(DrawCount), "Unexpected count type");

// 途中計算の優先精度浮動小数点の最低桁数(10進数)
using SizeOfDigits = unsigned int;
constexpr SizeOfDigits g_digitsOfFloat = 20;

// 途中計算の優先精度浮動小数点の最低桁数(2進数)
// 10進数の桁数 * { log_2(10) = ln(10)/ln(2) = 3.32193} 桁
constexpr SizeOfDigits GetBitsOfGmpFloat(void) {
    return static_cast<SizeOfDigits>(
        static_cast<double>(g_digitsOfFloat) * 3.32193 + 1.0);
}

constexpr CardNumber g_progressPerLoop = 100000; // 何種類ループを回したら進捗表示するか
constexpr int g_printedDoublePrecision = 15;     // doubleを表示するときの精度

// モンテカルロシミュレーション
class CardHolder {
    // ユニットテスト
    friend class CardUnitTest;
    friend class CombinationUnitTest;
public:
    CardHolder(CardNumber kindOfCards, RandomNumberSeed randomSeed);
    virtual ~CardHolder(void);
    virtual void Draw(CardNumber cardsPerDraw);
    virtual DrawCount Complete(CardNumber cardsPerDraw);
private:
    using HasCardArray = std::vector<bool>;
    virtual void fetchCardSet(CardNumber cardsPerDraw, CardNumberArray& cardArray);
    const CardNumber kindOfCards_;         // カードの全種類数
    DrawCount    drawCount_;               // セットを買った回数
    CardNumber   numberOfCollectedCards_;  // 何種類揃ったか
    RandomNumberGen  randGen_;   // 乱数発生器
    HasCardArray hasCardArray_;  // 各カードが揃ったかどうか
};

// 組み合わせ計算の基底クラス兼Factory
class CardProbabilityCalculator {
    // ユニットテスト
    friend class CardUnitTest;
    friend class CombinationUnitTest;
public:
    static CardProbabilityCalculator* CreateInstance(CardNumber kindOfCards,
                                                     CardNumber cardsPerDraw, Precision precision);
    CardProbabilityCalculator(CardNumber kindOfCards, CardNumber cardsPerDraw);
    virtual ~CardProbabilityCalculator(void);
    virtual ReportedCount Calculate(std::ostream* os);
    virtual const char* GetName(void) const = 0;
protected:
    virtual ReportedCount calculate(void) = 0;
    virtual void setProbabilityArray(CardNumber numberOfWantedCards) = 0;
    virtual void updateCount(CardNumber numberOfWantedCards) = 0;
    // 敢えて派生クラスから見えるようにしている
    const CardNumber kindOfCards_;    // カードの全種類数
    const CardNumber cardsPerDraw_;   // 1セットの種類数
    const CardNumber countArraySize_;        // カードの全種類についての配列の要素数
    const CardNumber probabilityArraySize_;  // カードのセットについての配列の要素数
};

// 浮動小数の型を指定した組み合わせ計算
template <typename FloatType>
class TypedCardProbabilityCaluculator : public CardProbabilityCalculator {
    // ユニットテスト
    friend class CardUnitTest;
    friend class CombinationUnitTest;
public:
    TypedCardProbabilityCaluculator(CardNumber kindOfCards, CardNumber cardsPerDraw);
    virtual ~TypedCardProbabilityCaluculator(void);
    virtual const char* GetName(void) const;
private:
    static const char* name_;  // モード名
    virtual ReportedCount calculate(void);
    virtual void setProbabilityArray(CardNumber numberOfWantedCards);
    virtual void updateCount(CardNumber numberOfWantedCards);
    virtual FloatType getCombinationSize(CardNumber setSize, CardNumber subSetSize);
    std::vector<FloatType> countFloatArray_;        // カードの残り種類数について回数の期待値
    std::vector<FloatType> probabilityFloatArray_;  // カードのセットについて遷移確率
    std::vector<FloatType> combinationFloatArray_;  // カードのセットについて遷移組み合わせ数
};

// ユニットテスト
class CardUnitTest {
public:
    CardUnitTest(void);
    virtual ~CardUnitTest(void);
    virtual int Run(void);
private:
    virtual int testCardHolder(void);
    virtual int testCardHolderComplete(void);
    virtual int testCardHolderFecth(void);
};

class CombinationUnitTest {
public:
    CombinationUnitTest(void);
    virtual ~CombinationUnitTest(void);
    virtual int Run(void);
private:
    virtual int testCalculator(void);
    virtual int testCalculatorGetCombinationSize(void);
};

/*
Local Variables:
mode: c++
coding: utf-8-unix
tab-width: nil
c-file-style: "stroustrup"
End:
*/
