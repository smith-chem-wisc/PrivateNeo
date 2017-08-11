using NeoInternal;
using NUnit.Framework;

namespace NeoTesting
{
    [TestFixture]
    public class NeoTests
    {
        [Test]
        public void test()
        {
            Assert.IsTrue(AlternativeSequences.ionsUsed.Count > 0);
        }
    }
}
