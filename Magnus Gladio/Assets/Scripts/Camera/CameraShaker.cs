using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CameraShaker : MonoBehaviour
{
    public IEnumerator Shake(float duration, float intensity)
    {
        Vector3 originalPosition = transform.position;

        float elapsedTime = 0;

        while (elapsedTime < duration)
        {
            float x = Random.Range(-1, 1) * intensity;
            float y = Random.Range(-1, 1) * intensity;

            transform.localPosition = new Vector3(x, y, originalPosition.z);

            elapsedTime += Time.unscaledDeltaTime;

            yield return null;
        }


        transform.position = originalPosition;

    }
}
