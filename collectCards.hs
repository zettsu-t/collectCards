-- 詳細はcollectCards.cpp を参照

import Data.List
import System.Environment (getArgs)

-- n枚からk枚を選ぶ組み合わせの総数
choose :: Int -> Int -> Double
choose n k | n <= 0 = 0
           | k < 0 = 0
           | k > n = 0
           | k == 0 = 1
           | n < (k + k) = choose n (n - k)
           | otherwise = (choose n (k - 1)) * (fromIntegral (n + 1 - k)) / (fromIntegral k)

-- 全n種類、残りr種類で、一度にd枚引くときに、0..d枚それぞれについて、
-- 引いたものの組み合わせ * 引けなかった(すでに持っていた)ものの組み合わせ
chooseVec n r d = map f [0..d]
  where f i = (choose r i) * (choose (n - r) (d - i))

-- 長さを1にする
toUnitVec ls = map (/ (sum ls)) ls

-- 全n種類、残りr種類で、一度にd枚引くときに 0..d枚増える確率
-- d枚に重複がないので、残り1種類ならd/nになる。なぜなら、
-- p(i=0) : 1 choose 0 * n-1 choose d = n-1 choose d
-- p(i=1) : 1 choose 1 * n-1 choose d-1 = p(i=0) * d / (n-d)
-- p(i>1) : 1 choose (>1) = 0
-- よって p(i=1) / p(i=0..d) = d / n
probVec n r d = toUnitVec $ chooseVec n r d

-- 全n種類、一度にd枚引くときに、残りr種類のときの回数の期待値を求める
nextProbVec n r d ls = [p] ++ (take d ls)
  where ps = probVec n r d
        qs = toUnitVec $ (drop 1 ps)
        p = (1.0 / (1.0 - ps!!0)) + sum (zipWith (*) qs ls)

-- 全n種類、一度にd枚引くときに、最初の一手(残りn-d種類)までの期待値を返す
lastProbVec n r d ls | d > (n - r) = ls
                     | otherwise = lastProbVec n (r + 1) d $ nextProbVec n r d ls

-- 最初の一手は別に数える
expectedCount n d = (f n d)!!0 + 1.0
  where f n d = lastProbVec n 1 d $ replicate d 0.0

-- 手数と価格を返す
solve n d price = (show count) ++ " (" ++ (show (count * price)) ++ ")"
  where count = expectedCount n d

-- コマンドラインから値を指定できるが、省略時の値も用意する
paramSet [] = (10, 1, 1.0)
paramSet (n:d:price:xs) = (read n :: Int, read d :: Int , read price :: Double)
paramSet (n:d:xs) = (read n :: Int, read d :: Int, 1.0)
paramSet (n:xs) = (read n :: Int, 1, 1.0)

main = do
  args <- getArgs
  let (n,d,price) = paramSet args
  putStrLn $ solve n d price
