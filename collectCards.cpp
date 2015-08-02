/*
 * カードが合計N種類あり、1セット(購入単位)には重複なしのカードがD枚
 * 入っているとする。このときN種類のカードが揃うために、必要なセット数の
 * 期待値を求める。
 */

#include <algorithm>
#include <iomanip>
#include <memory>
#include <unordered_map>
#include <assert.h>
#include <math.h>
#include <memory.h>
#include <stdint.h>
#include <time.h>
#include <boost/lexical_cast.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include "collectCards.hpp"

/*
 * コマンドライン引数を指定して起動する。後にある引数から順に省略可。
 * collectCards 第1引数 第2引数 第3引数 第4引数 第5引数
 * 第1引数 : シミュレーション回数 >= 0
 * 第2引数 : カードの全種類数(N) > 0
 * 第3引数 : 一セットのカードの種類数(D) < カードの全種類数
 * 第4引数 : 一セットのカードの価格 > 0
 * 第5引数 : 演算精度のenum値(0以上で範囲内)
 *   0 : C++組み込みのfloat
 *   1 : C++組み込みのdouble
 *   2 : 任意精度浮動小数(精度はコンパイル時に固定する)
 */

// デフォルト値
namespace {
    constexpr TrialCount g_DefaultTrials = 1;           // シミュレーション回数
    constexpr CardNumber g_DefaultKindOfCards = 10;     // カードの全種類数
    constexpr CardNumber g_DefaultCardsPerDraw = 1;     // 一セットのカードの種類数
    constexpr ReportedPrice g_DefaultPricePerDraw = 2;  // 一セットのカードの単価
    constexpr Precision  g_DefaultPrecision = Precision::PRECISION_BUILTIN_DOUBLE;  // 演算精度
}

/* 自明な解を求めて検算に用いる
セットの枚数が1のときは、カードを1枚引いて新しいカードが手に入る確率は、
単に残った種類の割合の逆数である(全部でN種類, 残りn種類ならN/n)。
毎回カードを買う行為は互いに独立なので、期待値は加法性がある。
そのためセットがそろうまでに引く回数の期待値は、
  sum(N/n) for n=1..N
で求まる。

この確率は、open hashのエントリが埋まった時に、各エントリに要素が
いくつあるかの期待値がO(log N)、全エントリの要素数がO(N*log(N))で
あることと同じ(1/nを積分するとlog(n)になる)。
*/

/*
 * モンテカルロシミュレーション
 * サンプリング誤差があるが、アルゴリズムを深く考えずにプログラミングが
 * 簡単にできるので速報値を求めることができる。アルゴリズムの検算にも
 * 使える。
 */
CardHolder::CardHolder(CardNumber kindOfCards, RandomNumberSeed randomSeed) :
    kindOfCards_(kindOfCards), drawCount_(0), numberOfCollectedCards_(0) {
    randGen_.seed(randomSeed);

    // falseで初期化
    hasCardArray_.assign(kindOfCards, false);
    return;
}

CardHolder::~CardHolder(void) {
    return;
}

void CardHolder::Draw(CardNumber cardsPerDraw) {
    if (!cardsPerDraw) {
        return;
    }

    // 0で初期化
    CardNumberArray cardArray(cardsPerDraw, 0);
    fetchCardSet(cardsPerDraw, cardArray);

    // すでに引いた種類のカードは無視する
    for(CardNumber i=0; i < cardsPerDraw; ++i) {
        auto num = cardArray.at(i);
        if (num >= kindOfCards_) {
            std::cerr << "Card number out of range! " << num << std::endl;
            ::abort();
        }

        if (!hasCardArray_.at(num)) {
            ++numberOfCollectedCards_;
            hasCardArray_.at(num) = true;
        }
    }

    return;
}

DrawCount CardHolder::Complete(CardNumber cardsPerDraw) {
    // 揃うまで引き続ける
    do {
        Draw(cardsPerDraw);
        ++drawCount_;
    } while (numberOfCollectedCards_ < kindOfCards_);

    return drawCount_;
}

void CardHolder::fetchCardSet(CardNumber cardsPerDraw, CardNumberArray& cardArray) {
    // どの種類のカードを引いたか
    std::unordered_map<CardNumber, bool> hasDrawnMap;
    // 両端を含む一様乱数
    std::uniform_int_distribution<CardNumber> randomCard(0, kindOfCards_ - 1);

    // カードの種類が重複したら引き直す
    for(CardNumber i=0; i<cardsPerDraw; ++i) {
        CardNumber num = 0;
        do {
            num = randomCard(randGen_);
        } while(hasDrawnMap[num]);

        cardArray[i] = num;
        hasDrawnMap[num] = true;
    }

    return;
}

/*
動的計画法(dynamic programming)で解く。
残りN-n種類のときの、セット数の期待値はいくらか、というのを求める。

記号の導入:
x[C]y : x個からy個を取り出す組み合わせ = x! / (y! * (x-y)!)
任意のxについて以下の式が成り立つ(計算を高速化するためと、条件分岐を
減らすためにこう定義する)
+ x[C]0 = 1, x[C]x = 1
+ x[C]1 = x, x[C]x-1 = x
+ x[C]y = x[C]x-y for y when 0 <= y <= x つまりyとx-yを入れ替えられる
+ x[C]y = 0 for y when y < 0 or y > x
+ 0[C]y = 0 for any y

一般的に、残りn種類とする。
+ 残り1種類を引く組み合わせは、
 - 引いたものの組み合わせn[C]1
 - 引けなかった(すでに持っていた)ものの組み合わせN-n[C]D-1
+ 一般的にm種類を引く組み合わせをc(n,m)とおくと、
 - 引いたものの組み合わせn[C]m
 - 引けなかった(すでに持っていた)ものの組み合わせN-n[C]D-m

nを固定して、残り0..m種類をそれぞれ引く組み合わせc(n,i) for i=0..n が、
全組み合わせ sum(c(n,i)) for i=0..n に占める確率を求める。

カードがほとんど揃っている場合、n + i > Nになるが、そのときは、
N-n[C]i = 0になり、N種類を超えて集められないことと整合性が取れる。
同様に最初の1セットではN[C]0 = 1, N-D[C]0 = 1 が成り立ち、必ずD種類
揃うことと整合性が取れる。

これを動的計画法で解く。たとえばN=100, D=2の場合、
残り100種類 <- 残り99種類
  : 残りの1種類を引くまでの回数の期待値 = 50 (100/2)
  セットに重複がないことに注意する。
  重複があるときは、(1/100)^2になり、まったく違う答えが出る。
100 <<- 98
    <- 99 <- 98
  : 残りの1または2種類を引くまでの回数の期待値
  + 99種類揃っていて残りの1種類を引くまでの期待値(もう求まっている) *
    (1回で残りのうち1種類を引く確率 / 1回で残りのうち1または2種類を引く確率)
99 <<- 97
   <- 98 <- 97
  : 残りの1または2種類を引くまでの回数の期待値
  + 99種類揃っていて残りの1種類を引くまでの期待値(もう求まっている) *
    (1回で残りのうち2種類を引く確率 / 1回で残りのうち1または2種類を引く確率)
  + 98種類揃っていて残りの1種類を引くまでの期待値(もう求まっている) *
    (1回で残りのうち1種類を引く確率 / 1回で残りのうち1または2種類を引く確率)
以下同様

残りn種類の時の回数の期待値 =
  n種類から一種類も増えない確率の逆数(=回数の期待値)
  + sum{(残りi種類増える確率/残り1..D種類増える確率)
    * n-i回のときの回数の期待値} for i=1..D
になる。これを残りN-D(最初の1セットを買ったとき)まで繰り返す。
*/

namespace {
    // Floatの桁数(10進数)
    using MultiPrecisionFloat = boost::multiprecision::number
        <boost::multiprecision::cpp_dec_float<g_digitsOfFloat>>;
    // Floatを十分格納できる桁数にしておく
    using MultiPrecisionInt = boost::multiprecision::checked_int512_t;
    using CPCalculatorFloat = TypedCardProbabilityCalculator<float>;
    using CPCalculatorDouble = TypedCardProbabilityCalculator<double>;
    using CPCalculatorMultiPrecision = TypedCardProbabilityCalculator<MultiPrecisionFloat>;
}

// Factory
CardProbabilityCalculator* CardProbabilityCalculator::CreateInstance(
    CardNumber kindOfCards, CardNumber cardsPerDraw, Precision precision) {
    switch(precision) {
    case Precision::PRECISION_BUILTIN_FLOAT:
        return new CPCalculatorFloat(kindOfCards, cardsPerDraw);
    case Precision::PRECISION_BUILTIN_DOUBLE:
        return new CPCalculatorDouble(kindOfCards, cardsPerDraw);
    case Precision::PRECISION_VARIABLE_PRECISION_FLOAT:
        return new CPCalculatorMultiPrecision(kindOfCards, cardsPerDraw);
    default:
        std::cerr << "Unknown precision enum" << static_cast<int>(precision) << std::endl;
        ::abort();
        break;
    }

    return nullptr;
}

CardProbabilityCalculator::CardProbabilityCalculator(CardNumber kindOfCards, CardNumber cardsPerDraw)
    : kindOfCards_(kindOfCards), cardsPerDraw_(cardsPerDraw),
      countArraySize_(kindOfCards + 1), probabilityArraySize_(cardsPerDraw + 1) {
    // 0..セットのカード数について確保するので1要素多く確保しておく
    return;
}

CardProbabilityCalculator::~CardProbabilityCalculator(void) {
    return;
}

ReportedCount CardProbabilityCalculator::Calculate(std::ostream* os) {
    // 残り1種類から、最初の1セットまで降る
    for(CardNumber i=1; i<=(kindOfCards_ - cardsPerDraw_); ++i) {
        if ((i % g_progressPerLoop) == 0) {
            if (os) {
                (*os) << "*" << std::flush;
            }
        }
        setProbabilityArray(i);
        updateCount(i);
    }

    if (os) {
        (*os) << std::endl;
    }

    // 最初の1手は別に数える
    constexpr ReportedCount firstDraw = 1;

    // template method
    return calculate() + firstDraw;
}

/* 型ごとに異なる定義 */
template <typename T> const char* TypedCardProbabilityCalculator<T>::name_ = "unknown type";
template <> const char* TypedCardProbabilityCalculator<float>::name_ = "builtin float";
template <> const char* TypedCardProbabilityCalculator<double>::name_ = "builtin double";
template <> const char* TypedCardProbabilityCalculator<MultiPrecisionFloat>::name_ = "variable precision float";

template <typename FloatType>
ReportedCount TypedCardProbabilityCalculator<FloatType>::calculate(void) {
    return countFloatArray_.at(0);
}

template <>
ReportedCount TypedCardProbabilityCalculator<MultiPrecisionFloat>::calculate(void) {
    return countFloatArray_.at(0).convert_to<ReportedCount>();
}

template <typename FloatType>
TypedCardProbabilityCalculator<FloatType>::TypedCardProbabilityCalculator(
    CardNumber kindOfCards, CardNumber cardsPerDraw) :
    CardProbabilityCalculator(kindOfCards, cardsPerDraw) {
    countFloatArray_.assign(countArraySize_, 0);
    probabilityFloatArray_.assign(probabilityArraySize_, 0);
    combinationFloatArray_.assign(probabilityArraySize_, 0);
}

template <typename FloatType>
TypedCardProbabilityCalculator<FloatType>::~TypedCardProbabilityCalculator(void) {
    return;
}

template <typename FloatType>
const char* TypedCardProbabilityCalculator<FloatType>::GetName(void) const {
    return name_;
}

template <typename FloatType>
void TypedCardProbabilityCalculator<FloatType>::setProbabilityArray(CardNumber numberOfWantedCards) {
    // 残り種類が少ないときは組み合わせが0になる
    for(CardNumber i=0; i<probabilityArraySize_; ++i) {
        probabilityFloatArray_.at(i) = 0;
        combinationFloatArray_.at(i) = 0;
    }

    // 残り種類が多すぎたら、1セットの種類より多くは集まらない
    const CardNumber numberOfCards = std::min(numberOfWantedCards, cardsPerDraw_);
    FloatType totalCombinationSize = 0;

    // 何種類増えるかの組み合わせの数を求める
    CombinationSizeSet comboSizeSet1;
    CombinationSizeSet comboSizeSet2;
    getCombinationSizeSet(numberOfWantedCards, numberOfCards, comboSizeSet1);
    getCombinationSizeSet(kindOfCards_ - numberOfWantedCards, cardsPerDraw_, comboSizeSet2);

    for(CardNumber i=0; i<=numberOfCards; ++i) {
        combinationFloatArray_.at(i) = comboSizeSet1.at(i) * comboSizeSet2.at(cardsPerDraw_ - i);
        totalCombinationSize += combinationFloatArray_.at(i);
    }

    // 何種類増えるかの確率分布を求める
    for(CardNumber i=0; i<=numberOfCards; ++i) {
        probabilityFloatArray_.at(i) = combinationFloatArray_.at(i) / totalCombinationSize;
    }

    return;
}

template <typename FloatType>
void TypedCardProbabilityCalculator<FloatType>::updateCount(CardNumber numberOfWantedCards) {
    // 種類が増える確率 = 1.0 - 1種類も増えない確率:probabilityArray_[0]
    // 種類が増えるまでの回数の期待値 = 確率の逆数
    const FloatType probabilityToForward = 1.0 - probabilityFloatArray_.at(0);
    FloatType totalExpectedCount = 1.0 / probabilityToForward;

    // 回数の期待値を、前回の計算から1種類分ずらす
    for(CardNumber i=cardsPerDraw_; i>=1; --i) {
        countFloatArray_.at(i) = countFloatArray_.at(i-1);
    }

    // ある種類数が増える確率に、ある種類増えたあとの回数期待値を掛ける
    for(CardNumber i=1; i<=cardsPerDraw_; ++i) {
        totalExpectedCount += probabilityFloatArray_.at(i) * countFloatArray_.at(i) / probabilityToForward;
    }

    // 回数の期待値をずらしてあるので先頭に入れる
    countFloatArray_.at(0) = totalExpectedCount;
    return;
}

// x[C]0, x[C]1, ... を一度に求めたほうが速い
// https://gist.github.com/keitaoouchi/2381408
// にある、大内慶太氏のScala実装からアイデアを頂きました。
template <typename FloatType>
void TypedCardProbabilityCalculator<FloatType>::getCombinationSizeSet(
    CardNumber setSize, CardNumber maxSubSetSize, CombinationSizeSet& comboSizeSet) {

    comboSizeSet.assign(maxSubSetSize+1, 0);

    for(CardNumber i=0; i<=maxSubSetSize; ++i) {
        // x[C]y = 0 for y when y < 0 or y > x , 0[C]y = 0 for any y
        if ((i > setSize) || (setSize == 0)) {
            comboSizeSet.at(i) = 0;
            continue;
        }

        // x[C]0 = 1, x[C]x = 1
        if ((i == 0) || (i == setSize)) {
            comboSizeSet.at(i) = 1;
            continue;
        }

        // x[C]1 = x, x[C]x-1 = x
        if ((i == 1) || ((i + 1) == setSize)) {
            comboSizeSet.at(i) = setSize;
            continue;
        }

        // x[C]i = x[C](i-1) * (x+1-i) / i
        comboSizeSet.at(i) = comboSizeSet.at(i-1);
        comboSizeSet.at(i) *= (setSize + 1 - i);
        comboSizeSet.at(i) /= i;
    }

    return;
}

namespace {
    // 指定回数だけシミュレーション
    ReportedCount simulate(TrialCount trials, CardNumber kindOfCards,
                           CardNumber cardsPerDraw, ReportedPrice pricePerDraw) {
        if (trials == 0) {
            return 0;
        }

        ReportedCount averangeCount = 0;
        DrawCount totalCount = 0;
        for(TrialCount i=1; i<=trials; ++i) {
            RandomNumberSeed randomSeed = static_cast<decltype(randomSeed)>(::time(nullptr) + i);
            CardHolder cardHolder(kindOfCards, randomSeed);
            totalCount += cardHolder.Complete(cardsPerDraw);
        }

        averangeCount = static_cast<decltype(averangeCount)>(totalCount) /
            static_cast<decltype(averangeCount)>(trials);
        ReportedCount averangeCards = static_cast<decltype(averangeCount)>(cardsPerDraw) * averangeCount;
        ReportedPrice averangePrice = averangeCount * pricePerDraw;

        std::cout << "simulated result (set)   : ";
        std::cout << std::setprecision(g_printedDoublePrecision) << averangeCount << "\n";
        std::cout << "simulated result (card)  : ";
        std::cout << std::setprecision(g_printedDoublePrecision) << averangeCards << "\n";
        std::cout << "simulated result (price) : ";
        std::cout << std::setprecision(g_printedDoublePrecision) << averangePrice << "\n";
        std::cout << "---\n" << std::flush;

        return averangeCount;
    }

    // 計算で求める
    ReportedCount calculate(CardNumber kindOfCards, CardNumber cardsPerDraw,
                            ReportedPrice pricePerDraw, Precision argPrecision) {
        ReportedCount calculatedCount = 0;
        // enumを++できない
        for(Precision precision=Precision::PRECISION_MIN; precision<Precision::PRECISION_MAX_PLUS_ONE;
            precision = static_cast<decltype(precision)>(1+static_cast<int>(precision))) {
            if (precision > argPrecision) {
                break;
            }

            std::unique_ptr<CardProbabilityCalculator> pCalc
                (CardProbabilityCalculator::CreateInstance(kindOfCards, cardsPerDraw, precision));
            std::cout << "Precision : " << pCalc->GetName() << "\n";
            calculatedCount = pCalc->Calculate(&std::cerr);
            const ReportedCount numberOfCards = static_cast<decltype(numberOfCards)>(cardsPerDraw) *
                calculatedCount;
            const ReportedPrice calculatedPrice = calculatedCount * pricePerDraw;

            std::cout << "calculated result (set)   : ";
            std::cout << std::setprecision(g_printedDoublePrecision) << calculatedCount << "\n";
            std::cout << "calculated result (card)  : ";
            std::cout << std::setprecision(g_printedDoublePrecision) << numberOfCards << "\n";
            std::cout << "calculated result (price) : ";
            std::cout << std::setprecision(g_printedDoublePrecision) << calculatedPrice << "\n";
            std::cout << "---\n\n" << std::flush;
        }

        return calculatedCount;
    }

    void exec(int argc, const char* const argv[]) {
        TrialCount trials = g_DefaultTrials;
        CardNumber kindOfCards = g_DefaultKindOfCards;
        CardNumber cardsPerDraw = g_DefaultCardsPerDraw;
        ReportedPrice pricePerDraw = g_DefaultPricePerDraw;
        Precision precision = g_DefaultPrecision;

        // 型変換に失敗したら即終了する
        if (argc > 1) {
            trials = boost::lexical_cast<decltype(trials)>(argv[1]);
        }

        if (argc > 2) {
            kindOfCards = boost::lexical_cast<decltype(kindOfCards)>(argv[2]);
        }

        if (argc > 3) {
            cardsPerDraw = boost::lexical_cast<decltype(cardsPerDraw)>(argv[3]);
        }

        if (argc > 4) {
            const auto argPrice = boost::lexical_cast<decltype(pricePerDraw)>(argv[4]);
            if (argPrice >= 1) {
                pricePerDraw = argPrice;
            }
        }

        if (argc > 5) {
            precision = static_cast<decltype(precision)>(boost::lexical_cast<int>(argv[5]));
            // とりあえず作ってみる。失敗したら終了。
            std::unique_ptr<CardProbabilityCalculator> pCalc
                (CardProbabilityCalculator::CreateInstance(3, 1, precision));
        }

        if ((kindOfCards <= 0) || (kindOfCards <= cardsPerDraw)) {
            std::cerr << "invalid argument\n";
            return;
        }

        std::cout << "\nkind of cards " << kindOfCards << ", cards per draw " << cardsPerDraw;
        std::cout << ", price to draw once " << pricePerDraw << "\n" << std::flush;

        const auto averangeCount = simulate(trials, kindOfCards, cardsPerDraw, pricePerDraw);
        const auto calculatedCount = calculate(kindOfCards, cardsPerDraw, pricePerDraw, precision);

        // サンプリング誤差を求める
        if (trials > 0) {
            std::cout << "sampling error rate [%] : ";
            std::cout << 100.0 * (averangeCount - calculatedCount) / calculatedCount << "\n";
            std::cout << std::flush;
        }

        return;
    }
}

CardUnitTest::CardUnitTest(void) {
    return;
}

CardUnitTest::~CardUnitTest(void) {
    return;
}

int CardUnitTest::Run(void) {
    return testCardHolder() + testCardHolderComplete() + testCardHolderFecth();
}

int CardUnitTest::testCardHolder(void) {
    constexpr CardNumber kindOfCards = 8;
    CardHolder cardHolder(kindOfCards, 1);

    assert(cardHolder.kindOfCards_ == kindOfCards);
    assert(cardHolder.drawCount_ == 0);
    assert(cardHolder.numberOfCollectedCards_ == 0);
    assert(cardHolder.hasCardArray_.size() == kindOfCards);
    assert(cardHolder.hasCardArray_.at(0) == false);
    assert(cardHolder.hasCardArray_.at(kindOfCards - 1) == false);

    // 本当のテストは再現可能な乱数生成モックに置き換える
    cardHolder.Draw(1);
    int count = 0;
    for(CardNumber i=0; i<kindOfCards; ++i) {
        if (cardHolder.hasCardArray_.at(i) == true) {
            ++count;
        }
    }
    assert(count == 1);

    return 0;
}

int CardUnitTest::testCardHolderComplete(void) {
    for(CardNumber kindOfCards = 2; kindOfCards < 10; ++kindOfCards) {
        CardHolder cardHolder(kindOfCards, 0xffffffff);
        for(CardNumber cardsPerDraw = 1; cardsPerDraw < kindOfCards - 1; ++cardsPerDraw) {
            auto count = cardHolder.Complete(cardsPerDraw);
            assert(cardHolder.numberOfCollectedCards_ == kindOfCards);
            assert(count >= (kindOfCards / cardsPerDraw));
        }
    }

    return 0;
}

int CardUnitTest::testCardHolderFecth(void) {
    for(CardNumber kindOfCards = 2; kindOfCards < 10; ++kindOfCards) {
        CardHolder cardHolder(kindOfCards, 0xffff);
        CardNumberArray cardArray(kindOfCards, 0);
        cardHolder.fetchCardSet(kindOfCards, cardArray);

        CardNumber totalNumber = 0;
        std::unordered_map<CardNumber, bool> hasDrawnMap;
        for(auto cardNumber : cardArray) {
            assert(cardNumber < kindOfCards);
            assert(hasDrawnMap[cardNumber] == false);
            totalNumber += cardNumber;
            hasDrawnMap[cardNumber] = true;
        }
        assert(totalNumber == (kindOfCards * (kindOfCards - 1) / 2));
    }

    return 0;
}

CombinationUnitTest::CombinationUnitTest(void) {
    return;
}

CombinationUnitTest::~CombinationUnitTest(void) {
    return;
}

int CombinationUnitTest::Run(void) {
    return testCalculator() + testCalculatorGetCombinationSize();
}

int CombinationUnitTest::testCalculator(void) {
    for(CardNumber kindOfCards = 2; kindOfCards < 10; ++kindOfCards) {
        CPCalculatorMultiPrecision calc(kindOfCards, 1);
        auto actual = calc.Calculate(0);
        decltype(actual) expected = 0;
        for(CardNumber i = 1; i <= kindOfCards; ++i) {
            expected += static_cast<decltype(expected)>(kindOfCards) /
                static_cast<decltype(expected)>(i);
        }
        assert(::fabs(actual - expected) < 0.001);
    }

    return 0;
}

int CombinationUnitTest::testCalculatorGetCombinationSize(void) {
    CPCalculatorMultiPrecision calc(100, 1);

//  constexpr CardNumber MaxSetSize = 11;  // uint32の範囲ではこれより大きいとオーバーフローする
    constexpr CardNumber MaxSetSize = 50;  // MultiPrecisionIntとMultiPrecisionFloatの限界

    for(CardNumber setSize = 0; setSize <= MaxSetSize; ++setSize) {
        const CardNumber upperSubSetSize = setSize + 2;
        decltype(calc)::CombinationSizeSet actualComboSet;
        calc.getCombinationSizeSet(setSize, upperSubSetSize, actualComboSet);

        MultiPrecisionInt expected;
        for(CardNumber subSetSize = 0; subSetSize <= upperSubSetSize; ++subSetSize) {
            if ((subSetSize > setSize) || (setSize == 0)) {
                expected = 0;
            } else if ((subSetSize == 0) || (subSetSize == setSize)) {
                expected = 1;
            } else {
                expected *= (setSize + 1 - subSetSize);
                expected /= subSetSize;
            }

            MultiPrecisionInt actual {MultiPrecisionFloat(
                    boost::multiprecision::floor(actualComboSet.at(subSetSize))).str()};
            assert(expected == actual);
        }
    }

    return 0;
}

int main(int argc, char* argv[]) {
    CardUnitTest test;
    assert(test.Run() == 0);

    CombinationUnitTest testCombo;
    assert(testCombo.Run() == 0);

    exec(argc, argv);
    return 0;
}

/*
Local Variables:
mode: c++
coding: utf-8-unix
tab-width: nil
c-file-style: "stroustrup"
End:
*/
